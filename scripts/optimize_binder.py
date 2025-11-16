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
from typing import Optional, Dict, List
import json
import numpy as np
import logging
import time

# BioPython for structure parsing
try:
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.SeqUtils import seq1
except ImportError:
    print("Warning: BioPython not installed. Sequence extraction will fail.")
    print("Install with: pip install biopython")

# Import optimization and output modules
sys.path.insert(0, str(Path(__file__).parent.parent))
from src.optimization.pareto import ParetoOptimizer
from src.output.json_generator import JSONOutputGenerator
from src.visualization.plots import SelectivityVisualizer
from src.selectivity.scorer import SelectivityScorer


class BinderOptimizationPipeline:
    """Orchestrates the complete binder optimization workflow."""

    def __init__(self, config_path: Path):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.project_name = self.config['project']['name']
        self.project_type = self.config['project']['type']
        self.output_dir = Path(self.config['output']['directory'])
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.start_time = datetime.now()

        # Setup logging
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = self.output_dir / f"pipeline_{self.timestamp}.log"
        self._setup_logging()

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

    def _setup_logging(self):
        """Setup logging to both file and console."""
        # Create formatters
        file_formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        console_formatter = logging.Formatter('%(message)s')

        # File handler
        file_handler = logging.FileHandler(self.log_file, mode='w')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)

        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(console_formatter)

        # Root logger
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

        self.logger = logger
        self.logger.info(f"Pipeline log file: {self.log_file}")

    def _log_stage_time(self, stage_name: str, start_time: float):
        """Log the time taken for a stage."""
        elapsed = time.time() - start_time
        mins = elapsed / 60
        self.logger.info(f"  ⏱ {stage_name} completed in {elapsed:.1f}s ({mins:.1f} min)")

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
        stage_start = time.time()
        print("\n[Stage 3/7] Running BoltzGen design generation...")

        params = self.config['parameters']
        num_designs = params.get('num_designs', 100)
        budget = params.get('budget', 50)

        workbench_dir = self.output_dir / "workbench"

        # Check if BoltzGen already completed
        metrics_file = workbench_dir / "intermediate_designs_inverse_folded" / "aggregate_metrics_analyze.csv"
        if metrics_file.exists():
            print(f"  ⚠ BoltzGen results already exist at: {workbench_dir}")
            print(f"  ✓ Skipping BoltzGen generation (using existing results)")
            self._log_stage_time("Stage 3: BoltzGen (skipped)", stage_start)
            return workbench_dir

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
            self._log_stage_time("Stage 3: BoltzGen", stage_start)
            return workbench_dir
        except subprocess.CalledProcessError as e:
            print(f"\nERROR: BoltzGen failed with exit code {e.returncode}")
            sys.exit(1)

    def _extract_sequences_from_structure(self, structure_path: Path, chain_ids: list = None) -> Dict[str, str]:
        """
        Extract amino acid sequences from PDB/CIF structure files.

        Args:
            structure_path: Path to PDB or CIF file
            chain_ids: List of chain IDs to extract (if None, extracts all)

        Returns:
            Dict mapping chain_id -> amino acid sequence
        """
        sequences = {}

        try:
            # Choose parser based on file extension
            if structure_path.suffix.lower() == '.cif':
                parser = MMCIFParser(QUIET=True)
            else:
                parser = PDBParser(QUIET=True)

            structure = parser.get_structure('struct', str(structure_path))

            # Extract sequences from each chain
            for model in structure:
                for chain in model:
                    chain_id = chain.get_id()

                    # Skip if specific chains requested and this isn't one
                    if chain_ids is not None and chain_id not in chain_ids:
                        continue

                    # Build sequence from residues
                    seq_list = []
                    for residue in chain:
                        # Skip non-standard residues (water, ligands, etc.)
                        if residue.id[0] == ' ':  # Standard residue
                            try:
                                # Convert 3-letter code to 1-letter
                                aa = seq1(residue.get_resname())
                                seq_list.append(aa)
                            except KeyError:
                                # Unknown residue, skip
                                pass

                    if seq_list:
                        sequences[chain_id] = ''.join(seq_list)

            return sequences

        except Exception as e:
            print(f"    ⚠ Error extracting sequences from {structure_path}: {e}")
            return {}

    def _run_offtarget_predictions(self, workbench_dir: Path) -> Optional[Path]:
        """
        Stage 3.5: Run BoltzGen predictions for off-target binding.

        For each top N candidate × each off-target:
        - Create YAML for [Design + Off-target] complex
        - Run BoltzGen folding prediction
        - Save to workbench/offtarget_predictions/

        Configuration:
        - num_offtarget_predictions: How many top candidates to score against off-targets
          * If not specified: uses num_candidates (smart pre-filtering, recommended)
          * If "all": scores all num_designs (exhaustive, slower)
          * If number: scores that many top candidates
        """
        stage_start = time.time()
        print("\n[Stage 3.5/7] Running off-target predictions...")

        # Check if off-targets are specified
        off_targets = self.config['targets'].get('negative', [])
        if not off_targets:
            print("  ⚠ No off-targets specified, skipping off-target predictions")
            return None

        params = self.config['parameters']

        # Determine how many candidates to score against off-targets
        # Option B (default): Smart pre-filtering - score only top candidates
        # Option A: Score all designs - set num_offtarget_predictions = "all" or num_designs
        num_offtarget_predictions = params.get('num_offtarget_predictions', None)

        if num_offtarget_predictions == "all":
            # Option A: Score all designs
            num_candidates = params.get('num_designs', 100)
            print(f"  Mode: Score ALL designs against off-targets (exhaustive)")
        elif num_offtarget_predictions is not None:
            # Custom number specified
            num_candidates = int(num_offtarget_predictions)
            print(f"  Mode: Score top {num_candidates} candidates against off-targets (custom)")
        else:
            # Option B (default): Smart pre-filtering using num_candidates
            num_candidates = params.get('num_candidates', 30)
            print(f"  Mode: Smart pre-filtering - top {num_candidates} candidates only (recommended)")

        # Load primary metrics to get top candidates
        # BoltzGen saves metrics in intermediate_designs_inverse_folded/ subdirectory
        metrics_file = workbench_dir / "intermediate_designs_inverse_folded" / "aggregate_metrics_analyze.csv"
        if not metrics_file.exists():
            print(f"  ⚠ Metrics file not found: {metrics_file}")
            print("  Skipping off-target predictions")
            return None

        df = pd.read_csv(metrics_file)
        if 'design_to_target_iptm' not in df.columns:
            print("  ⚠ Missing design_to_target_iptm column, skipping off-target predictions")
            return None

        # Get top N candidates
        top_candidates = df.sort_values('design_to_target_iptm', ascending=False).head(num_candidates)
        print(f"  Selected top {len(top_candidates)} candidates for off-target prediction")

        # Create output directory
        offtarget_predictions_dir = workbench_dir / "offtarget_predictions"
        offtarget_predictions_dir.mkdir(exist_ok=True)

        # Find design structures directory
        design_dirs = [
            workbench_dir / "designs",
            workbench_dir / "inverse_fold_cifs",
            workbench_dir / "design_cifs",
            workbench_dir / "intermediate_designs",
            workbench_dir / "intermediate_designs_inverse_folded"
        ]

        structure_dir = None
        for d in design_dirs:
            if d.exists():
                structure_dir = d
                break

        if structure_dir is None:
            print(f"  ⚠ No design structures found in {workbench_dir}")
            return None

        print(f"  Using structures from: {structure_dir}")

        # Process each off-target
        for off_target in off_targets:
            off_target_name = off_target['name']
            off_target_pdb = Path(off_target['pdb'])

            print(f"\n  [{off_targets.index(off_target)+1}/{len(off_targets)}] Processing off-target: {off_target_name}")

            # Create directory for this off-target
            offtarget_dir = offtarget_predictions_dir / off_target_name
            offtarget_dir.mkdir(exist_ok=True)

            # Check if predictions already exist for this off-target
            # Count how many prediction directories exist
            existing_predictions = 0
            for idx, row in top_candidates.iterrows():
                design_id = row.get('id', row.get('design_id', f'design_{idx}'))
                prediction_dir = offtarget_dir / design_id
                # Check for Boltz-2 output structure
                yaml_stem = f"{design_id}_vs_{off_target_name}"
                boltz_results = prediction_dir / f"boltz_results_{yaml_stem}"
                confidence_file = boltz_results / "predictions" / yaml_stem / f"confidence_{yaml_stem}_model_0.json"
                if confidence_file.exists():
                    existing_predictions += 1

            if existing_predictions == len(top_candidates):
                print(f"    ✓ All {existing_predictions} predictions already exist, skipping {off_target_name}")
                continue
            elif existing_predictions > 0:
                print(f"    ⚠ Found {existing_predictions}/{len(top_candidates)} existing predictions, will complete remaining...")

            # Get off-target chain ID
            offtarget_chain = off_target.get('chain', 'A')

            # Process each candidate
            for idx, row in top_candidates.iterrows():
                design_id = row.get('id', row.get('design_id', f'design_{idx}'))

                # Check if prediction already exists for this design
                prediction_output = offtarget_dir / design_id
                yaml_stem = f"{design_id}_vs_{off_target_name}"
                boltz_results = prediction_output / f"boltz_results_{yaml_stem}"
                confidence_file = boltz_results / "predictions" / yaml_stem / f"confidence_{yaml_stem}_model_0.json"

                if confidence_file.exists():
                    # Prediction already exists, skip
                    continue

                # Find structure file
                structure_file = None
                for ext in ['.cif', '.pdb']:
                    candidate_file = structure_dir / f"{design_id}{ext}"
                    if candidate_file.exists():
                        structure_file = candidate_file
                        break

                if structure_file is None:
                    print(f"    ⚠ Structure not found for {design_id}")
                    continue

                # Extract sequences from structures
                # BUG FIX: Only extract BINDER chains, not the primary target chain
                # BoltzGen outputs: first N chains = designed binder, last chain = primary target
                # where N = len(scaffold.chains)
                num_scaffold_chains = len(self.config['scaffold']['chains'])

                all_design_sequences = self._extract_sequences_from_structure(structure_file)
                if not all_design_sequences:
                    print(f"    ⚠ Failed to extract sequences for {design_id}")
                    continue

                # Only take the first N chains (the designed binder chains)
                all_chain_ids = sorted(all_design_sequences.keys())  # Sort to ensure order
                binder_chain_ids = all_chain_ids[:num_scaffold_chains]
                design_sequences = {chain_id: all_design_sequences[chain_id] for chain_id in binder_chain_ids}

                offtarget_sequences = self._extract_sequences_from_structure(off_target_pdb, [offtarget_chain])

                if not design_sequences or not offtarget_sequences:
                    print(f"    ⚠ Failed to extract sequences for {design_id}")
                    continue

                # Get design chain IDs
                design_chains = list(design_sequences.keys())

                # Assign unique chain IDs for the complex
                # Binder chains: keep original IDs (A, B, C, etc.)
                # Off-target chain: use next available ID
                all_chain_ids = set(design_chains)
                offtarget_new_id = offtarget_chain
                if offtarget_chain in all_chain_ids:
                    # Find next available chain ID
                    for candidate_id in 'DEFGHIJKLMNOPQRSTUVWXYZ':
                        if candidate_id not in all_chain_ids:
                            offtarget_new_id = candidate_id
                            break

                # Create Boltz-2 compatible YAML with sequences
                sequences_list = []

                # Add binder sequences - IMPORTANT: create separate entry for each chain!
                # This ensures Boltz-2 treats them as different proteins (heterodimer)
                # NOT as a homotrimer
                for chain_id in design_chains:
                    sequences_list.append({
                        'protein': {
                            'id': chain_id,
                            'sequence': design_sequences[chain_id]
                        }
                    })

                # Add off-target sequence
                sequences_list.append({
                    'protein': {
                        'id': offtarget_new_id,
                        'sequence': offtarget_sequences[offtarget_chain]
                    }
                })

                yaml_content = {
                    'version': 1,
                    'sequences': sequences_list
                }

                # Save YAML
                yaml_file = offtarget_dir / f"{design_id}_vs_{off_target_name}.yaml"
                with open(yaml_file, 'w') as f:
                    yaml.dump(yaml_content, f)

                # Run Boltz-2 prediction directly
                prediction_output.mkdir(exist_ok=True)

                # Use boltz predict with MSA server
                cmd_predict = [
                    'boltz', 'predict',
                    str(yaml_file),
                    '--out_dir', str(prediction_output),
                    '--devices', '1',  # Use GPU 1
                    '--accelerator', 'gpu',
                    '--recycling_steps', '3',
                    '--sampling_steps', '200',
                    '--model', 'boltz2',  # Use Boltz-2 model
                    '--use_msa_server',  # Auto-generate MSAs
                    '--override'
                ]

                try:
                    # Run Boltz-2 prediction
                    result = subprocess.run(cmd_predict, capture_output=True, timeout=600, text=True)

                    if result.returncode != 0:
                        print(f"    ⚠ Boltz-2 failed for {design_id}: {result.stderr[:200]}")
                        continue

                    # Extract ipTM and PAE from the prediction outputs
                    metrics_dict = self._extract_iptm_from_boltz2(prediction_output, yaml_file.stem)

                    if metrics_dict is not None:
                        # Save metrics for this prediction
                        self._save_prediction_metrics(prediction_output, design_id, off_target_name, metrics_dict)
                        iptm = metrics_dict.get('iptm', 0.0)
                        pae = metrics_dict.get('pae')
                        pae_str = f", PAE = {pae:.2f}" if pae is not None else ""
                        print(f"      ✓ {design_id}: ipTM = {iptm:.3f}{pae_str}")
                    else:
                        print(f"    ⚠ Could not extract metrics for {design_id}")

                except subprocess.TimeoutExpired:
                    print(f"    ⚠ Timeout for {design_id} (>10min)")
                    continue
                except Exception as e:
                    print(f"    ⚠ Error for {design_id}: {e}")
                    continue

            print(f"    ✓ Completed {off_target_name} predictions")

        print(f"\n  ✓ Off-target predictions saved to: {offtarget_predictions_dir}")
        self._log_stage_time("Stage 3.5: Off-target predictions", stage_start)
        return offtarget_predictions_dir

    def _get_chain_ids(self, structure_file: Path) -> List[str]:
        """Extract chain IDs from a structure file."""
        import gemmi

        st = gemmi.read_structure(str(structure_file))
        chain_ids = []

        for model in st:
            for chain in model:
                if chain.name not in chain_ids:
                    chain_ids.append(chain.name)

        return chain_ids

    def _extract_iptm_from_boltz2(self, prediction_dir: Path, yaml_stem: str) -> Optional[Dict[str, float]]:
        """
        Extract ipTM and PAE scores from Boltz-2 prediction outputs.

        Boltz-2 output structure (actual):
        prediction_dir/
          boltz_results_{yaml_stem}/
            predictions/
              {yaml_stem}/
                confidence_{yaml_stem}_model_0.json  <- Contains ipTM
                pae_{yaml_stem}_model_0.npz          <- Contains PAE matrix

        Returns:
            Dict with 'iptm' and 'pae' keys, or None if not found
        """
        import glob

        # Primary path: look for confidence JSON file in actual Boltz-2 structure
        boltz_results_dir = prediction_dir / f"boltz_results_{yaml_stem}"
        confidence_file = boltz_results_dir / "predictions" / yaml_stem / f"confidence_{yaml_stem}_model_0.json"

        if confidence_file.exists():
            try:
                with open(confidence_file, 'r') as f:
                    data = json.load(f)
                    result = {}

                    # Extract ipTM
                    if 'iptm' in data:
                        result['iptm'] = float(data['iptm'])

                    # Try to load PAE from separate .npz file
                    pae_file = boltz_results_dir / "predictions" / yaml_stem / f"pae_{yaml_stem}_model_0.npz"
                    if pae_file.exists():
                        try:
                            pae_npz = np.load(pae_file)
                            # PAE matrix is typically stored under 'pae' key
                            if 'pae' in pae_npz:
                                pae_matrix = pae_npz['pae']
                                result['pae'] = float(np.mean(pae_matrix))
                            else:
                                # Try first array in npz
                                pae_matrix = pae_npz[pae_npz.files[0]]
                                result['pae'] = float(np.mean(pae_matrix))
                        except Exception as e:
                            print(f"      ⚠ Could not load PAE from npz: {e}")

                    if 'iptm' in result:
                        return result
            except Exception as e:
                print(f"      ⚠ Error reading confidence file: {e}")

        # Fallback: search for any confidence JSON file
        boltz_pred_dir = prediction_dir / yaml_stem / "predictions"
        if boltz_pred_dir.exists():
            json_pattern = str(boltz_pred_dir / "*_confidence*.json")
            json_files = glob.glob(json_pattern)

            for json_file in json_files:
                try:
                    with open(json_file, 'r') as f:
                        data = json.load(f)
                        result = {}

                        if 'iptm' in data:
                            result['iptm'] = float(data['iptm'])

                        if 'pae' in data:
                            pae_data = data['pae']
                            if isinstance(pae_data, (list, np.ndarray)):
                                result['pae'] = float(np.mean(pae_data))
                            else:
                                result['pae'] = float(pae_data)

                        if 'iptm' in result:
                            return result
                except Exception:
                    continue

        # Last resort: try alternative output structures
        alt_paths = [
            prediction_dir / "predictions" / yaml_stem / "predictions",
            prediction_dir / "predictions",
        ]

        for alt_dir in alt_paths:
            if alt_dir.exists():
                json_pattern = str(alt_dir / "*_confidence*.json")
                json_files = glob.glob(json_pattern)
                for json_file in json_files:
                    try:
                        with open(json_file, 'r') as f:
                            data = json.load(f)
                            result = {}

                            if 'iptm' in data:
                                result['iptm'] = float(data['iptm'])

                            if 'pae' in data:
                                pae_data = data['pae']
                                if isinstance(pae_data, (list, np.ndarray)):
                                    result['pae'] = float(np.mean(pae_data))
                                else:
                                    result['pae'] = float(pae_data)

                            if 'iptm' in result:
                                return result
                    except Exception:
                        continue

        return None

    def _save_prediction_metrics(self, prediction_dir: Path, design_id: str, offtarget_name: str, metrics_dict: Dict[str, float]):
        """Save prediction metrics to CSV for the scorer to read."""
        import pandas as pd

        # Create a simple metrics file that the scorer can read
        metrics_file = prediction_dir / "prediction_metrics.csv"

        iptm_score = metrics_dict.get('iptm', 0.0)
        pae_score = metrics_dict.get('pae', 20.0)

        metrics = {
            'design_id': [f"{design_id}_vs_{offtarget_name}"],
            'design_to_target_iptm': [iptm_score],
            'iptm': [iptm_score],
            'min_design_to_target_pae': [pae_score],
            'pae': [pae_score]
        }

        df = pd.DataFrame(metrics)
        df.to_csv(metrics_file, index=False)

        pae_str = f", PAE = {pae_score:.2f}" if 'pae' in metrics_dict else ""
        print(f"      ✓ Extracted ipTM = {iptm_score:.3f}{pae_str}")

    def _run_selectivity_scoring(self, workbench_dir: Path, offtarget_predictions_dir: Optional[Path]):
        """Score designs for selectivity against off-targets."""
        print("\n[Stage 4/7] Running selectivity scoring...")

        params = self.config['parameters']
        num_candidates = params.get('num_candidates', 30)

        # Build primary and off-target lists
        primary_target = self.config['targets']['positive'][0]
        primary_pdb = Path(primary_target['pdb'])

        # Build off-targets list with (path, name) tuples using config names
        off_targets = []
        for target in self.config['targets'].get('negative', []):
            off_targets.append((Path(target['pdb']), target['name']))

        if not off_targets:
            print("  ⚠ No off-targets specified, skipping selectivity scoring")
            return None

        scorer_output = self.output_dir / "selectivity_results"
        scorer_output.mkdir(exist_ok=True, parents=True)

        print(f"  Scoring top {num_candidates} candidates...")
        print(f"  Primary target: {primary_target['name']}")
        print(f"  Off-targets: {', '.join([name for _, name in off_targets])}")
        if offtarget_predictions_dir:
            print(f"  Using pre-computed predictions from: {offtarget_predictions_dir}")

        try:
            # Initialize scorer with correct off-target names from config
            scorer = SelectivityScorer(
                primary_target=primary_pdb,
                off_targets=off_targets,
                offtarget_predictions_dir=offtarget_predictions_dir
            )

            # Run scoring
            results_df = scorer.score_variants(
                variants_dir=workbench_dir,
                output_dir=scorer_output,
                num_candidates=num_candidates
            )

            print("\n  ✓ Selectivity scoring complete")
            return scorer_output
        except Exception as e:
            print(f"\nERROR: Selectivity scoring failed: {e}")
            import traceback
            traceback.print_exc()
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
        total_time = time.time() - self.start_time.timestamp()
        total_mins = total_time / 60
        total_hours = total_mins / 60

        print("\n[Stage 7/7] Pipeline summary")
        print(f"  Project: {self.project_name}")
        print(f"  Type: {self.project_type}")
        print(f"  Output directory: {self.output_dir}")
        print(f"  Log file: {self.log_file}")
        print(f"  Total time: {total_time:.1f}s ({total_mins:.1f} min / {total_hours:.2f} hours)")

        self.logger.info(f"\n{'='*80}")
        self.logger.info(f"PIPELINE COMPLETED")
        self.logger.info(f"Total execution time: {total_time:.1f}s ({total_mins:.1f} min / {total_hours:.2f} hours)")
        self.logger.info(f"{'='*80}")

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

            # Stage 3.5: Run off-target predictions
            offtarget_predictions_dir = self._run_offtarget_predictions(workbench_dir)

            # Stage 4: Run selectivity scoring
            scorer_output = self._run_selectivity_scoring(workbench_dir, offtarget_predictions_dir)

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
