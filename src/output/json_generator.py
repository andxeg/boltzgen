"""
JSON Output Generator for Binder Optimization Results
======================================================

Generates structured JSON output for easy web UI integration.
"""

import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional


class JSONOutputGenerator:
    """Generate structured JSON output from optimization results."""

    def __init__(self, config: Dict):
        """
        Initialize JSON generator.

        Args:
            config: Original YAML configuration
        """
        self.config = config
        self.project_name = config['project']['name']
        self.project_type = config['project']['type']

    def generate_output(
        self,
        df: pd.DataFrame,
        pareto_results: Optional[Dict] = None,
        output_path: Optional[Path] = None
    ) -> Dict:
        """
        Generate complete JSON output.

        Args:
            df: DataFrame with all metrics and scores
            pareto_results: Results from ParetoOptimizer (if available)
            output_path: Path to save JSON file

        Returns:
            JSON-serializable dict
        """
        # Build JSON structure according to plan
        output = {
            "project": self._generate_project_metadata(),
            "summary": self._generate_summary(df, pareto_results),
            "ranked_designs": self._generate_ranked_designs(df, pareto_results),
            "pareto_analysis": self._generate_pareto_analysis(pareto_results) if pareto_results else None,
            "visualizations": self._list_visualizations(output_path.parent if output_path else None)
        }

        # Save to file
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_path, 'w') as f:
                json.dump(output, f, indent=2)

            print(f"  ✓ JSON output saved: {output_path}")

        return output

    def _generate_project_metadata(self) -> Dict:
        """Generate project metadata section."""
        positive_targets = [t['name'] for t in self.config['targets']['positive']]
        negative_targets = [t['name'] for t in self.config['targets'].get('negative', [])]

        return {
            "name": self.project_name,
            "type": self.project_type,
            "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "targets": {
                "positive": positive_targets,
                "negative": negative_targets
            }
        }

    def _generate_summary(
        self,
        df: pd.DataFrame,
        pareto_results: Optional[Dict] = None
    ) -> Dict:
        """Generate summary statistics."""
        summary = {
            "total_designs": int(len(df)),
            "scored_designs": int(len(df)),
        }

        # Add Pareto stats if available
        if pareto_results:
            summary["pareto_optimal"] = int(pareto_results['num_pareto'])
        else:
            summary["pareto_optimal"] = 0

        # Add mean metrics
        if 'design_to_target_iptm' in df.columns:
            summary["mean_affinity"] = float(df['design_to_target_iptm'].mean())

        if 'selectivity_composite' in df.columns:
            summary["mean_selectivity"] = float(df['selectivity_composite'].mean())

        if 'min_design_to_target_pae' in df.columns:
            summary["mean_pae"] = float(df['min_design_to_target_pae'].mean())

        return summary

    def _generate_ranked_designs(
        self,
        df: pd.DataFrame,
        pareto_results: Optional[Dict] = None,
        top_n: int = 30
    ) -> List[Dict]:
        """Generate ranked designs list."""
        # Sort by selectivity score or pareto rank
        if 'pareto_rank' in df.columns:
            df_sorted = df.sort_values(['pareto_rank', 'selectivity_composite'], ascending=[True, False])
        elif 'selectivity_composite' in df.columns:
            df_sorted = df.sort_values('selectivity_composite', ascending=False)
        else:
            df_sorted = df.sort_values('design_to_target_iptm', ascending=False)

        # Take top N
        df_top = df_sorted.head(top_n)

        designs = []
        for rank, (idx, row) in enumerate(df_top.iterrows(), 1):
            design = {
                "rank": rank,
                "design_id": str(row['design_id']),
                "sequence": str(row.get('sequence', '')),
                "metrics": self._extract_metrics(row),
                "composite_score": float(row.get('selectivity_composite', 0.0)),
                "pareto_optimal": bool(row.get('pareto_optimal', False)),
                "pareto_rank": int(row.get('pareto_rank', -1)) if 'pareto_rank' in row else None
            }

            # Add structure file path if available
            if 'structure_file' in row:
                design["structure_file"] = str(row['structure_file'])

            designs.append(design)

        return designs

    def _extract_metrics(self, row: pd.Series) -> Dict:
        """Extract and organize metrics for a single design."""
        metrics = {}

        # Affinity metrics
        affinity = {}
        if 'design_to_target_iptm' in row:
            affinity['primary_iptm'] = float(row['design_to_target_iptm'])
        if 'min_design_to_target_pae' in row:
            affinity['primary_pae'] = float(row['min_design_to_target_pae'])

        if affinity:
            metrics['affinity'] = affinity

        # Selectivity metrics
        selectivity = {}
        off_target_cols = [col for col in row.index if col.endswith('_iptm') and col.startswith('offtarget_')]

        if off_target_cols:
            # Extract off-target names and scores
            for col in off_target_cols:
                target_name = col.replace('offtarget_', '').replace('_iptm', '')
                selectivity[f"{target_name}_iptm"] = float(row[col])

        if 'max_offtarget_iptm' in row:
            selectivity['max_offtarget_iptm'] = float(row['max_offtarget_iptm'])

        if 'selectivity_iptm' in row:
            selectivity['selectivity_score'] = float(row['selectivity_iptm'])

        if selectivity:
            metrics['selectivity'] = selectivity

        # Biophysical properties
        properties = {}
        if 'delta_sasa_refolded' in row:
            properties['sasa'] = float(row['delta_sasa_refolded'])
        if 'solubility' in row:
            properties['solubility'] = float(row['solubility'])
        if 'stability' in row:
            properties['stability'] = float(row['stability'])
        if 'design_ptm' in row:
            properties['ptm'] = float(row['design_ptm'])

        if properties:
            metrics['properties'] = properties

        return metrics

    def _generate_pareto_analysis(self, pareto_results: Dict) -> Dict:
        """Generate Pareto analysis section."""
        if not pareto_results:
            return None

        analysis = pareto_results.get('analysis', {})

        # Get Pareto frontier design IDs
        frontier_df = pareto_results['df'][pareto_results['df']['pareto_optimal']]
        frontier_ranks = frontier_df.sort_values('selectivity_composite', ascending=False).index.tolist()
        frontier_ranks = [int(r) + 1 for r in frontier_ranks]  # Convert to 1-indexed ranks

        pareto_output = {
            "frontier_designs": frontier_ranks[:10],  # Top 10 Pareto designs
            "num_pareto": int(pareto_results['num_pareto']),
            "trade_offs": {}
        }

        # Add trade-off descriptions
        if 'trade_offs' in analysis:
            pareto_output['trade_offs'] = analysis['trade_offs']

        # Add correlation info
        if 'correlations' in analysis:
            pareto_output['correlations'] = {
                k: float(v) for k, v in analysis['correlations'].items()
            }

        return pareto_output

    def _list_visualizations(self, output_dir: Optional[Path]) -> List[str]:
        """List available visualization files."""
        if not output_dir:
            return []

        plots_dir = output_dir / "plots"
        if not plots_dir.exists():
            return []

        vis_files = []
        for plot_file in plots_dir.glob("*.png"):
            vis_files.append(f"plots/{plot_file.name}")

        return sorted(vis_files)


def generate_json_output(
    results_csv: Path,
    config: Dict,
    pareto_csv: Optional[Path] = None,
    output_path: Optional[Path] = None
) -> Dict:
    """
    Convenience function to generate JSON output.

    Args:
        results_csv: Path to selectivity results CSV
        config: Original YAML configuration
        pareto_csv: Path to Pareto-optimized CSV (optional)
        output_path: Path to save JSON file

    Returns:
        JSON output dict
    """
    # Load results
    df = pd.read_csv(results_csv)

    # Load Pareto results if available
    pareto_results = None
    if pareto_csv and pareto_csv.exists():
        pareto_df = pd.read_csv(pareto_csv)

        # Merge Pareto columns
        if 'pareto_optimal' in pareto_df.columns:
            df = pareto_df

        pareto_results = {
            'df': df,
            'num_pareto': int(df['pareto_optimal'].sum()) if 'pareto_optimal' in df.columns else 0,
            'analysis': {}  # Can be populated from separate analysis file
        }

    # Generate output
    generator = JSONOutputGenerator(config)
    return generator.generate_output(df, pareto_results, output_path)


if __name__ == '__main__':
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description='Generate JSON output from optimization results')
    parser.add_argument('--results', type=Path, required=True, help='Selectivity results CSV')
    parser.add_argument('--config', type=Path, required=True, help='Original YAML config')
    parser.add_argument('--pareto', type=Path, help='Pareto-optimized CSV (optional)')
    parser.add_argument('--output', type=Path, required=True, help='Output JSON path')

    args = parser.parse_args()

    # Load config
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # Generate JSON
    json_output = generate_json_output(
        args.results,
        config,
        args.pareto,
        args.output
    )

    print(f"✓ JSON output generated: {args.output}")
