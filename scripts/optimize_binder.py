#!/usr/bin/env python3
"""
Portable Binder Optimization Pipeline
=====================================

Main entry point for the fine-grained binder design platform.
Supports both existing binder optimization and de novo design with selectivity guidance.

Usage:
    python scripts/optimize_binder.py config/cxcr4_optimization.yaml
"""

import argparse
import sys
import yaml
import subprocess
from pathlib import Path
from datetime import datetime
import shutil
import pandas as pd
from typing import Optional, Dict

# Import optimization and output modules
sys.path.insert(0, str(Path(__file__).parent.parent))
from src.optimization.pareto import ParetoOptimizer
from src.output.json_generator import JSONOutputGenerator
from src.visualization.plots import SelectivityVisualizer


class BinderOptimizationPipeline:
    """Orchestrates the complete binder optimization workflow."""

    def __init__(self, config_path: Path):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.project_name = self.config['project']['name']
        self.project_type = self.config['project']['type']
        self.output_dir = Path(self.config['output']['directory'])
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    def _load_config(self) -> dict:
        """Load and validate YAML configuration."""
        try:
            with open(self.config_path, 'r') as f:
                config = yaml.safe_load(f)

            # Validate required fields
            required_keys = ['project', 'targets', 'parameters', 'output']
            for key in required_keys:
                if key not in config:
                    raise ValueError(f"Missing required config section: {key}")

            # Check project type
            if config['project']['type'] not in ['existing_binder', 'de_novo']:
                raise ValueError(f"Invalid project type: {config['project']['type']}")

            # Validate targets
            if 'positive' not in config['targets']:
                raise ValueError("Missing positive targets")

            print(f"✓ Configuration loaded: {self.config_path}")
            return config

        except Exception as e:
            print(f"ERROR: Failed to load config: {e}")
            sys.exit(1)

    def _validate_input_files(self):
        """Verify all input files exist."""
        print("\n[Stage 1/7] Validating input files...")

        missing_files = []

        # Check positive target PDBs
        for target in self.config['targets']['positive']:
            pdb_path = Path(target['pdb'])
            if not pdb_path.exists():
                missing_files.append(str(pdb_path))

        # Check negative target PDBs
        for target in self.config['targets'].get('negative', []):
            pdb_path = Path(target['pdb'])
            if not pdb_path.exists():
                missing_files.append(str(pdb_path))

        # Check scaffold PDB (for existing binder)
        if self.project_type == 'existing_binder':
            if 'scaffold' not in self.config:
                raise ValueError("Missing scaffold configuration for existing_binder mode")
            scaffold_pdb = Path(self.config['scaffold']['pdb'])
            if not scaffold_pdb.exists():
                missing_files.append(str(scaffold_pdb))

        if missing_files:
            print("ERROR: Missing input files:")
            for f in missing_files:
                print(f"  - {f}")
            sys.exit(1)

        print("  ✓ All input files found")

    def _generate_boltzgen_config(self) -> Path:
        """Generate BoltzGen YAML configuration."""
        print("\n[Stage 2/7] Generating BoltzGen configuration...")

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Generate BoltzGen config based on project type
        if self.project_type == 'existing_binder':
            boltzgen_config = self._generate_inverse_folding_config()
        else:
            boltzgen_config = self._generate_denovo_config()

        # Save config
        config_path = self.output_dir / "boltzgen_config.yaml"
        with open(config_path, 'w') as f:
            yaml.dump(boltzgen_config, f, default_flow_style=False)

        print(f"  ✓ BoltzGen config saved: {config_path}")
        return config_path

    def _generate_inverse_folding_config(self) -> dict:
        """Generate config for existing binder optimization (inverse folding)."""
        scaffold = self.config['scaffold']
        primary_target = self.config['targets']['positive'][0]

        # Convert paths to absolute
        scaffold_path = Path(scaffold['pdb']).resolve()
        target_path = Path(primary_target['pdb']).resolve()

        # Build design regions for BoltzGen entities format
        design_spec = []
        for region in scaffold.get('design_regions', []):
            # Convert "26-35,50-65,95-102" to "26..35,50..65,95..102"
            res_index = region['residues'].replace('-', '..')
            design_spec.append({
                'chain': {
                    'id': region['chain'],
                    'res_index': res_index
                }
            })

        # Build include spec for scaffold chains
        include_spec = []
        for chain in scaffold.get('chains', []):
            include_spec.append({
                'chain': {
                    'id': chain
                }
            })

        # BoltzGen entities format
        config = {
            'entities': [
                # Scaffold structure with design regions
                {
                    'file': {
                        'path': str(scaffold_path),
                        'include': include_spec,
                        'design': design_spec
                    }
                },
                # Primary target
                {
                    'file': {
                        'path': str(target_path),
                        'include': [
                            {
                                'chain': {
                                    'id': primary_target.get('chain', 'A')
                                }
                            }
                        ]
                    }
                }
            ]
        }

        return config

    def _generate_denovo_config(self) -> dict:
        """Generate config for de novo binder design."""
        primary_target = self.config['targets']['positive'][0]

        # Convert path to absolute
        target_path = Path(primary_target['pdb']).resolve()

        # For de novo mode, we only need the target - BoltzGen will design the binder
        # BoltzGen entities format
        config = {
            'entities': [
                # Primary target only - binder will be designed de novo
                {
                    'file': {
                        'path': str(target_path),
                        'include': [
                            {
                                'chain': {
                                    'id': primary_target.get('chain', 'R')
                                }
                            }
                        ]
                    }
                }
            ]
        }

        return config

    def _run_boltzgen(self, config_path: Path):
        """Execute BoltzGen design generation."""
        print("\n[Stage 3/7] Running BoltzGen design generation...")

        params = self.config['parameters']
        num_designs = params.get('num_designs', 100)
        budget = params.get('budget', 50)

        workbench_dir = self.output_dir / "workbench"

        # Build BoltzGen command
        cmd = [
            'boltzgen', 'run',
            str(config_path),
            '--output', str(workbench_dir),
            '--protocol', 'nanobody-anything',
            '--num_designs', str(num_designs),
            '--budget', str(budget),
            '--devices', '1'
        ]

        print(f"  Command: {' '.join(cmd)}")
        print(f"  Generating {num_designs} designs with budget {budget}...")
        print(f"  This may take some time...\n")

        # Run BoltzGen
        try:
            result = subprocess.run(cmd, check=True, capture_output=False)
            print("\n  ✓ BoltzGen generation complete")
            return workbench_dir
        except subprocess.CalledProcessError as e:
            print(f"\nERROR: BoltzGen failed with exit code {e.returncode}")
            sys.exit(1)

    def _run_selectivity_scoring(self, workbench_dir: Path):
        """Score designs for selectivity against off-targets."""
        print("\n[Stage 4/7] Running selectivity scoring...")

        params = self.config['parameters']
        num_candidates = params.get('num_candidates', 30)

        # Build primary and off-target lists
        primary_target = self.config['targets']['positive'][0]
        primary_pdb = primary_target['pdb']

        offtarget_pdbs = []
        for target in self.config['targets'].get('negative', []):
            offtarget_pdbs.append(target['pdb'])

        if not offtarget_pdbs:
            print("  ⚠ No off-targets specified, skipping selectivity scoring")
            return None

        # Build scorer command
        scorer_output = self.output_dir / "selectivity_results"

        cmd = [
            'python', 'src/selectivity/scorer.py',
            '--variants', str(workbench_dir),
            '--primary', primary_pdb,
            '--off_targets'
        ] + offtarget_pdbs + [
            '--output', str(scorer_output),
            '--num_candidates', str(num_candidates)
        ]

        print(f"  Scoring top {num_candidates} candidates...")
        print(f"  Primary target: {primary_target['name']}")
        print(f"  Off-targets: {', '.join([t['name'] for t in self.config['targets']['negative']])}")

        try:
            result = subprocess.run(cmd, check=True, cwd=Path.cwd())
            print("\n  ✓ Selectivity scoring complete")
            return scorer_output
        except subprocess.CalledProcessError as e:
            print(f"\nERROR: Selectivity scoring failed with exit code {e.returncode}")
            sys.exit(1)

    def _run_pareto_optimization(self, scorer_output: Path) -> Optional[Path]:
        """Run Pareto multi-objective optimization."""
        scoring_config = self.config.get('scoring', {})

        if not scoring_config.get('multi_objective', False):
            print("\n[Stage 5/7] Pareto optimization: Skipped (multi_objective=false)")
            return None

        print("\n[Stage 5/7] Running Pareto multi-objective optimization...")

        csv_file = scorer_output / "aggregate_metrics_selectivity.csv"
        if not csv_file.exists():
            print("  ⚠ Selectivity CSV not found, skipping Pareto optimization")
            return None

        # Load scored data
        df = pd.read_csv(csv_file)

        # Get weights from config
        weights = scoring_config.get('weights', {
            'affinity': 0.6,
            'selectivity': 0.3,
            'properties': 0.1
        })

        # Initialize Pareto optimizer
        optimizer = ParetoOptimizer(
            objectives=['design_to_target_iptm', 'selectivity_composite'],
            maximize={'design_to_target_iptm': True, 'selectivity_composite': True}
        )

        # Run optimization
        pareto_results = optimizer.optimize(
            df,
            objectives=['design_to_target_iptm', 'selectivity_composite']
        )

        # Save Pareto results
        pareto_output = self.output_dir / "pareto_optimization"
        pareto_output.mkdir(exist_ok=True)

        # Extract Pareto-optimal designs from dataframe
        pareto_df = pareto_results['df'][pareto_results['df']['pareto_optimal']]

        pareto_df.to_csv(
            pareto_output / "pareto_frontier.csv",
            index=False
        )

        # Save trade-off analysis
        if 'analysis' in pareto_results:
            import json
            trade_off_file = pareto_output / "trade_off_analysis.json"
            with open(trade_off_file, 'w') as f:
                json.dump(pareto_results['analysis'], f, indent=2)

        print(f"  ✓ Pareto frontier: {len(pareto_df)} designs")
        print(f"  ✓ Results saved: {pareto_output}/")

        return pareto_results

    def _generate_outputs(self, scorer_output: Path, pareto_results: Optional[Dict]):
        """Generate final outputs in requested formats."""
        print("\n[Stage 6/7] Generating output files...")

        output_formats = self.config['output'].get('format', ['csv', 'json'])

        # Copy CSV results
        if 'csv' in output_formats and scorer_output:
            csv_src = scorer_output / "aggregate_metrics_selectivity.csv"
            if csv_src.exists():
                csv_dst = self.output_dir / "ranked_designs.csv"
                shutil.copy(csv_src, csv_dst)
                print(f"  ✓ CSV results: {csv_dst}")

        # Generate JSON output
        if 'json' in output_formats:
            print("  Generating JSON output...")
            csv_file = scorer_output / "aggregate_metrics_selectivity.csv"

            if csv_file.exists():
                # Load the dataframe (either with Pareto results or without)
                if pareto_results:
                    df = pareto_results['df']
                    pareto_dict = pareto_results
                else:
                    df = pd.read_csv(csv_file)
                    pareto_dict = None

                json_gen = JSONOutputGenerator(self.config)
                json_dst = self.output_dir / "results.json"

                json_output = json_gen.generate_output(
                    df=df,
                    pareto_results=pareto_dict,
                    output_path=json_dst
                )

                print(f"  ✓ JSON results: {json_dst}")
            else:
                print("  ⚠ CSV file not found, skipping JSON generation")

        # Generate visualization plots
        if 'plots' in output_formats:
            print("  Generating visualization plots...")
            csv_file = scorer_output / "aggregate_metrics_selectivity.csv"

            if csv_file.exists():
                plots_dst = self.output_dir / "plots"

                visualizer = SelectivityVisualizer(dpi=300)
                metrics_df = pd.read_csv(csv_file)

                visualizer.create_comprehensive_dashboard(
                    metrics_df,
                    plots_dst,
                    title_prefix=self.project_name
                )

                visualizer.create_summary_table(metrics_df, plots_dst)

                print(f"  ✓ Visualizations: {plots_dst}/")
            else:
                print("  ⚠ CSV file not found, skipping visualization generation")

        print(f"\n✓ All outputs saved to: {self.output_dir}")

    def _print_final_summary(self):
        """Print final pipeline summary."""
        print("\n[Stage 7/7] Pipeline summary")
        print(f"  Project: {self.project_name}")
        print(f"  Type: {self.project_type}")
        print(f"  Output directory: {self.output_dir}")

    def run(self):
        """Execute the complete pipeline."""
        print("="*80)
        print(f"Fine-Grained Binder Design Pipeline")
        print("="*80)
        print(f"Project: {self.project_name}")
        print(f"Type: {self.project_type}")
        print(f"Output: {self.output_dir}")
        print("="*80)

        try:
            # Stage 1: Validate inputs
            self._validate_input_files()

            # Stage 2: Generate BoltzGen config
            boltzgen_config = self._generate_boltzgen_config()

            # Stage 3: Run BoltzGen
            workbench_dir = self._run_boltzgen(boltzgen_config)

            # Stage 4: Run selectivity scoring
            scorer_output = self._run_selectivity_scoring(workbench_dir)

            # Stage 5: Run Pareto optimization
            pareto_results = self._run_pareto_optimization(scorer_output) if scorer_output else None

            # Stage 6: Generate outputs
            self._generate_outputs(scorer_output, pareto_results)

            # Stage 7: Summary
            self._print_final_summary()

            print("\n" + "="*80)
            print("✓ PIPELINE COMPLETE!")
            print("="*80)
            print(f"\nResults available at: {self.output_dir}")
            print(f"  - ranked_designs.csv (Top candidates with selectivity scores)")
            print(f"  - results.json (Structured output for web UI)")
            print(f"  - pareto_optimization/ (Multi-objective analysis)")
            print(f"  - plots/ (Visualization)")
            print("="*80 + "\n")

        except Exception as e:
            print(f"\nERROR: Pipeline failed: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Optimize binder design with selectivity guidance'
    )
    parser.add_argument(
        'config',
        type=Path,
        help='Path to YAML configuration file'
    )

    args = parser.parse_args()

    if not args.config.exists():
        print(f"ERROR: Configuration file not found: {args.config}")
        sys.exit(1)

    # Run pipeline
    pipeline = BinderOptimizationPipeline(args.config)
    pipeline.run()


if __name__ == '__main__':
    main()
