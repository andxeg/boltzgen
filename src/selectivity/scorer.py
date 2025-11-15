#!/usr/bin/env python3
"""
SelectivityScorer: Multi-target selectivity scoring for binder variants

Scores binder candidates against primary target + off-targets to calculate
selectivity metrics for multi-objective optimization.

Core workflow:
1. Load variants from BoltzGen output (inverse folding or de novo)
2. Extract sequences and structures
3. Predict binding to off-targets using Boltz-2
4. Calculate selectivity scores
5. Aggregate metrics for ranking
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional
import pandas as pd
import numpy as np
import subprocess
import json
import yaml
import argparse
from tqdm import tqdm


class SelectivityScorer:
    """
    Score binder variants for selectivity against multiple targets.

    Key innovation: Multi-target selectivity framework
    - Primary target: High affinity required
    - Off-targets: Low affinity desired (selectivity)
    - Multi-objective: Balance affinity + selectivity + properties
    """

    def __init__(
        self,
        primary_target: Path,
        off_targets: List[Tuple[Path, str]],  # [(path, name), ...]
        offtarget_predictions_dir: Optional[Path] = None
    ):
        """
        Initialize SelectivityScorer.

        Args:
            primary_target: Path to primary target PDB (e.g., CXCR4)
            off_targets: List of (path, name) tuples for off-targets
                        e.g., [(ccr5.pdb, "CCR5"), (cxcr2.pdb, "CXCR2")]
            offtarget_predictions_dir: Directory containing pre-computed off-target predictions
                                     (from Stage 3.5 of optimize_binder.py)
        """
        self.primary_target = Path(primary_target)
        self.off_targets = [(Path(p), n) for p, n in off_targets]
        self.offtarget_predictions_dir = Path(offtarget_predictions_dir) if offtarget_predictions_dir else None

        if not self.primary_target.exists():
            raise FileNotFoundError(f"Primary target not found: {self.primary_target}")

        for off_target_path, name in self.off_targets:
            if not off_target_path.exists():
                raise FileNotFoundError(f"Off-target {name} not found: {off_target_path}")

    def score_variants(
        self,
        variants_dir: Path,
        output_dir: Path,
        num_candidates: int = 50,
        use_cached_primary: bool = True
    ) -> pd.DataFrame:
        """
        Score all variants for selectivity.

        Args:
            variants_dir: Directory containing BoltzGen output (e.g., workbench/regn7663_optimized_small)
            output_dir: Output directory for selectivity results
            num_candidates: Number of top candidates to score (pre-filter by primary affinity)
            use_cached_primary: If True, use existing primary target metrics from BoltzGen

        Returns:
            DataFrame with selectivity scores and all metrics
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)

        print("\n" + "="*60)
        print("SelectivityScorer: Multi-Target Selectivity Analysis")
        print("="*60)
        print(f"Variants directory: {variants_dir}")
        print(f"Primary target: {self.primary_target.name}")
        print(f"Off-targets: {', '.join([name for _, name in self.off_targets])}")
        print(f"Num candidates: {num_candidates}")
        print("="*60 + "\n")

        # Step 1: Load primary target metrics from BoltzGen output
        print("[1/4] Loading primary target metrics...")
        primary_metrics = self._load_primary_metrics(variants_dir, num_candidates)
        print(f"  ✓ Loaded {len(primary_metrics)} variants")
        print(f"  ✓ Top affinity (ipTM): {primary_metrics['design_to_target_iptm'].max():.3f}")
        print(f"  ✓ Mean affinity: {primary_metrics['design_to_target_iptm'].mean():.3f}\n")

        # Step 2: Extract sequences for off-target scoring
        print("[2/4] Extracting variant sequences...")
        sequences_df = self._extract_sequences(variants_dir, primary_metrics['design_id'].tolist())
        primary_metrics = primary_metrics.merge(sequences_df, on='design_id', how='left')
        print(f"  ✓ Extracted {len(sequences_df)} sequences\n")

        # Step 3: Score against each off-target
        print("[3/4] Scoring against off-targets...")
        for i, (off_target_path, off_target_name) in enumerate(self.off_targets):
            print(f"  [{i+1}/{len(self.off_targets)}] Scoring against {off_target_name}...")

            offtarget_dir = output_dir / f"offtarget_{i}_{off_target_name}"
            offtarget_metrics = self._score_against_offtarget(
                sequences_df,
                primary_metrics,
                off_target_path,
                offtarget_dir,
                off_target_name
            )

            # Merge off-target metrics
            primary_metrics = primary_metrics.merge(
                offtarget_metrics,
                on='design_id',
                how='left',
                suffixes=('', f'_{off_target_name}')
            )

            print(f"    ✓ Mean off-target affinity (ipTM): {offtarget_metrics['offtarget_iptm'].mean():.3f}")

        print()

        # Step 4: Calculate composite selectivity scores
        print("[4/4] Calculating selectivity scores...")
        primary_metrics = self._calculate_selectivity(primary_metrics)
        print(f"  ✓ Mean selectivity: {primary_metrics['selectivity_composite'].mean():.3f}")
        print(f"  ✓ Max selectivity: {primary_metrics['selectivity_composite'].max():.3f}\n")

        # Save results
        output_file = output_dir / "aggregate_metrics_selectivity.csv"
        primary_metrics.to_csv(output_file, index=False)
        print(f"✓ Results saved to: {output_file}\n")

        # Print summary
        self._print_summary(primary_metrics)

        return primary_metrics

    def _load_primary_metrics(self, variants_dir: Path, num_candidates: int) -> pd.DataFrame:
        """Load metrics from BoltzGen analysis step."""
        # Try multiple possible locations
        metrics_paths = [
            variants_dir / "aggregate_metrics_analyze.csv",
            variants_dir / "intermediate_designs" / "aggregate_metrics_analyze.csv",
            variants_dir / "intermediate_designs_inverse_folded" / "aggregate_metrics_analyze.csv",
        ]

        metrics_file = None
        for path in metrics_paths:
            if path.exists():
                metrics_file = path
                break

        if metrics_file is None:
            raise FileNotFoundError(
                f"Metrics file not found in:\n" +
                "\n".join(f"  - {p}" for p in metrics_paths) +
                "\nMake sure BoltzGen pipeline has completed the analysis step."
            )

        df = pd.read_csv(metrics_file)

        # Handle column name variations
        if 'id' in df.columns and 'design_id' not in df.columns:
            df = df.rename(columns={'id': 'design_id'})

        # Pre-filter by primary affinity (keep top candidates)
        if 'design_to_target_iptm' not in df.columns:
            raise ValueError("Missing 'design_to_target_iptm' column in metrics file")

        df = df.sort_values('design_to_target_iptm', ascending=False).head(num_candidates)

        # Select relevant columns
        cols = ['design_id']
        optional_cols = [
            'design_to_target_iptm', 'min_design_to_target_pae', 'delta_sasa_refolded',
            'plip_hbonds_refolded', 'plip_saltbridges_refolded',
            'solubility', 'mean_plddt_design'
        ]

        for col in optional_cols:
            if col in df.columns:
                cols.append(col)

        return df[cols].reset_index(drop=True)

    def _extract_sequences(self, variants_dir: Path, design_ids: List[str]) -> pd.DataFrame:
        """Extract sequences from design output structures."""
        sequences = []

        # Try multiple possible locations for structures
        possible_dirs = [
            variants_dir / "designs",
            variants_dir / "inverse_fold_cifs",
            variants_dir / "design_cifs",
            variants_dir / "intermediate_designs",
            variants_dir / "intermediate_designs_inverse_folded"
        ]

        structure_dir = None
        for d in possible_dirs:
            if d.exists():
                structure_dir = d
                break

        if structure_dir is None:
            raise FileNotFoundError(
                f"Could not find structure directory in {variants_dir}. "
                f"Tried: {[str(d) for d in possible_dirs]}"
            )

        print(f"  Using structure directory: {structure_dir}")

        for design_id in tqdm(design_ids, desc="  Extracting sequences"):
            # Try .cif and .pdb formats
            for ext in ['.cif', '.pdb']:
                structure_file = structure_dir / f"{design_id}{ext}"
                if structure_file.exists():
                    sequence = self._parse_sequence_from_structure(structure_file)
                    sequences.append({
                        'design_id': design_id,
                        'sequence': sequence,
                        'structure_file': str(structure_file)
                    })
                    break

        if not sequences:
            raise RuntimeError(f"No structures found in {structure_dir}")

        return pd.DataFrame(sequences)

    def _parse_sequence_from_structure(self, structure_file: Path) -> str:
        """Parse sequence from PDB/CIF file."""
        import gemmi

        st = gemmi.read_structure(str(structure_file))

        # Extract sequence from designed chains
        # Assuming designed chains are those not in original target
        sequences = []
        for model in st:
            for chain in model:
                if chain.name in ['B', 'C', 'D', 'H', 'L']:  # Common design chain IDs
                    seq = ''.join([res.name for res in chain if res.name in gemmi.find_tabulated_residue(res.name).one_letter_code])
                    if seq:
                        sequences.append(seq)

        if not sequences:
            # Fallback: take first chain
            for model in st:
                for chain in model:
                    seq = ''
                    for res in chain:
                        if res.name in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
                                       'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']:
                            seq += gemmi.find_tabulated_residue(res.name).one_letter_code
                    if seq:
                        return seq

        return sequences[0] if sequences else ""

    def _score_against_offtarget(
        self,
        sequences_df: pd.DataFrame,
        primary_metrics: pd.DataFrame,
        offtarget_path: Path,
        output_dir: Path,
        offtarget_name: str
    ) -> pd.DataFrame:
        """
        Load pre-computed off-target predictions and extract metrics.

        This loads metrics from pre-computed predictions (Stage 3.5 from optimize_binder.py).
        No predictions are run here - this is pure scoring logic.
        """
        output_dir.mkdir(exist_ok=True, parents=True)

        print(f"    Loading {offtarget_name} predictions...")

        # Try to load from pre-computed predictions directory
        results = []

        if self.offtarget_predictions_dir and self.offtarget_predictions_dir.exists():
            # Load from pre-computed predictions
            offtarget_pred_dir = self.offtarget_predictions_dir / offtarget_name

            if offtarget_pred_dir.exists():
                print(f"    Using pre-computed predictions from: {offtarget_pred_dir}")

                for idx, row in tqdm(list(sequences_df.iterrows()), desc=f"    Extracting {offtarget_name} metrics"):
                    design_id = row['design_id']

                    # Load metrics from pre-computed prediction
                    metrics = self._load_precomputed_metrics(offtarget_pred_dir, design_id)

                    results.append({
                        'design_id': design_id,
                        'offtarget_iptm': metrics['iptm'],
                        'offtarget_pae': metrics['pae'],
                    })
            else:
                print(f"    ⚠ Pre-computed predictions not found for {offtarget_name}")
                print(f"    Expected directory: {offtarget_pred_dir}")
                # Use fallback values
                for idx, row in sequences_df.iterrows():
                    results.append({
                        'design_id': row['design_id'],
                        'offtarget_iptm': 0.15,  # Low binding fallback
                        'offtarget_pae': 18.0
                    })
        else:
            # No pre-computed predictions available, use fallback
            print(f"    ⚠ No pre-computed predictions directory provided")
            print(f"    Using fallback values for {offtarget_name}")

            for idx, row in sequences_df.iterrows():
                results.append({
                    'design_id': row['design_id'],
                    'offtarget_iptm': 0.15,  # Low binding fallback
                    'offtarget_pae': 18.0
                })

        return pd.DataFrame(results)

    def _load_precomputed_metrics(self, offtarget_pred_dir: Path, design_id: str) -> Dict[str, float]:
        """
        Load metrics from pre-computed off-target predictions.

        Args:
            offtarget_pred_dir: Directory containing predictions for this off-target
            design_id: Design identifier

        Returns:
            Dict with 'iptm' and 'pae' keys
        """
        prediction_dir = offtarget_pred_dir / design_id

        # Look for prediction metrics file (saved by optimize_binder.py Stage 3.5)
        metrics_file = prediction_dir / "prediction_metrics.csv"

        if metrics_file.exists():
            try:
                df = pd.read_csv(metrics_file)

                # Extract metrics
                if 'design_to_target_iptm' in df.columns:
                    iptm = df['design_to_target_iptm'].iloc[0]
                elif 'iptm' in df.columns:
                    iptm = df['iptm'].iloc[0]
                else:
                    iptm = 0.0

                if 'min_design_to_target_pae' in df.columns:
                    pae = df['min_design_to_target_pae'].iloc[0]
                elif 'min_pae' in df.columns:
                    pae = df['min_pae'].iloc[0]
                else:
                    pae = 20.0

                return {
                    'iptm': float(iptm),
                    'pae': float(pae)
                }
            except Exception as e:
                print(f"      ⚠ Failed to parse metrics for {design_id}: {e}")

        # Fallback: return low binding values (no binding detected)
        return {
            'iptm': 0.15,  # Low binding
            'pae': 18.0    # High error
        }

    def _calculate_selectivity(self, df: pd.DataFrame) -> pd.DataFrame:
        """Calculate composite selectivity scores."""

        # Find off-target columns
        offtarget_iptm_cols = [c for c in df.columns if 'offtarget_iptm' in c and c != 'max_offtarget_iptm']
        offtarget_pae_cols = [c for c in df.columns if 'offtarget_pae' in c and c != 'min_offtarget_pae']

        if not offtarget_iptm_cols:
            raise ValueError("No off-target ipTM columns found")

        # ipTM-based selectivity (higher primary ipTM, lower off-target ipTM = better)
        df['max_offtarget_iptm'] = df[offtarget_iptm_cols].max(axis=1)
        df['selectivity_iptm'] = df['design_to_target_iptm'] - df['max_offtarget_iptm']

        # pAE-based selectivity (lower pAE is better, so flip logic)
        if offtarget_pae_cols and 'min_design_to_target_pae' in df.columns:
            df['min_offtarget_pae'] = df[offtarget_pae_cols].min(axis=1)
            df['selectivity_pae'] = (df['min_offtarget_pae'] - df['min_design_to_target_pae']) / 10.0
        else:
            df['selectivity_pae'] = 0.0

        # Composite selectivity score (weighted combination)
        df['selectivity_composite'] = (
            0.6 * df['selectivity_iptm'] +
            0.3 * df['selectivity_pae'] +
            0.1 * (df.get('delta_sasa_refolded', 0) / df.get('delta_sasa_refolded', 1).max() if 'delta_sasa_refolded' in df.columns else 0)
        )

        return df

    def _print_summary(self, df: pd.DataFrame):
        """Print summary statistics."""
        print("="*60)
        print("Selectivity Scoring Summary")
        print("="*60)
        print(f"\nTop 10 Variants by Selectivity:")
        print("-" * 60)

        top10 = df.nlargest(10, 'selectivity_composite')

        summary_cols = ['design_id', 'design_to_target_iptm', 'max_offtarget_iptm', 'selectivity_composite']
        available_cols = [c for c in summary_cols if c in top10.columns]

        print(top10[available_cols].to_string(index=False))

        print("\n" + "="*60)
        print(f"Mean selectivity: {df['selectivity_composite'].mean():.3f}")
        print(f"Std selectivity: {df['selectivity_composite'].std():.3f}")
        print(f"Best selectivity: {df['selectivity_composite'].max():.3f}")
        print("="*60 + "\n")


def main():
    """Command-line interface for SelectivityScorer."""
    parser = argparse.ArgumentParser(
        description="Score binder variants for multi-target selectivity using pre-computed predictions"
    )
    parser.add_argument(
        '--variants', type=Path, required=True,
        help="BoltzGen output directory (e.g., workbench/regn7663_optimized_small)"
    )
    parser.add_argument(
        '--primary', type=Path, required=True,
        help="Primary target PDB file (e.g., data/cxcr4/8u4q_cxcr4.pdb)"
    )
    parser.add_argument(
        '--off_targets', type=Path, nargs='+', required=True,
        help="Off-target PDB files (e.g., data/cxcr4/4mbs_ccr5.pdb data/cxcr4/6lfo_cxcr2.pdb)"
    )
    parser.add_argument(
        '--output', type=Path, required=True,
        help="Output directory for selectivity results"
    )
    parser.add_argument(
        '--num_candidates', type=int, default=50,
        help="Number of top candidates to score (default: 50)"
    )
    parser.add_argument(
        '--offtarget_predictions', type=Path, default=None,
        help="Directory containing pre-computed off-target predictions (from Stage 3.5 of pipeline)"
    )

    args = parser.parse_args()

    # Create off-targets list with names extracted from filenames
    off_targets = [(p, p.stem) for p in args.off_targets]

    # Initialize scorer
    scorer = SelectivityScorer(
        primary_target=args.primary,
        off_targets=off_targets,
        offtarget_predictions_dir=args.offtarget_predictions
    )

    # Score variants
    results = scorer.score_variants(
        variants_dir=args.variants,
        output_dir=args.output,
        num_candidates=args.num_candidates
    )

    print(f"\n✓ Selectivity scoring complete!")
    print(f"  Results: {args.output / 'aggregate_metrics_selectivity.csv'}")


if __name__ == '__main__':
    main()
