#!/usr/bin/env python3
"""
Simple scoring script that uses Boltzgen's pre-computed metrics.
No need to run IPSAE - Boltzgen already computed similar metrics during analysis.
"""

import os
import pandas as pd
from pathlib import Path

# --- CONFIGURATION ---
RESULT_DIRS = [
    "/lambda/nfs/Nipah-hackathon/andrei/boltzgen/protein_competition/results/20designs_5budget"
]

OUTPUT_CSV = "final_rankings_boltzgen.csv"


def load_boltzgen_metrics(result_dir):
    """
    Load metrics from Boltzgen's aggregate_metrics_analyze.csv file.
    This file contains pre-computed design_to_target_iptm, interaction_pae, etc.
    """
    result_path = Path(result_dir)

    # Check for the inverse folded results (highest quality)
    metrics_file = result_path / "intermediate_designs_inverse_folded" / "aggregate_metrics_analyze.csv"

    if not metrics_file.exists():
        # Fallback to regular intermediate designs
        metrics_file = result_path / "intermediate_designs" / "aggregate_metrics_analyze.csv"

    if not metrics_file.exists():
        print(f"Warning: Could not find metrics file in {result_dir}")
        return pd.DataFrame()

    print(f"Loading metrics from: {metrics_file}")
    df = pd.read_csv(metrics_file)

    # Add source experiment column
    df['source_experiment'] = result_path.name

    return df


def score_all_directories():
    """Score all result directories and combine results."""
    all_metrics = []

    for result_dir in RESULT_DIRS:
        if not os.path.exists(result_dir):
            print(f"Warning: Directory not found: {result_dir}")
            continue

        print(f"\n--- Processing {result_dir} ---")
        metrics_df = load_boltzgen_metrics(result_dir)

        if not metrics_df.empty:
            all_metrics.append(metrics_df)
            print(f"  ✓ Loaded {len(metrics_df)} candidates")

    if not all_metrics:
        print("\nNo metrics found!")
        return None

    # Combine all metrics
    combined_df = pd.concat(all_metrics, ignore_index=True)

    # Select key columns for ranking
    ranking_columns = [
        'id',
        'source_experiment',
        'designed_sequence',
        'num_design',
        # Key affinity metrics (higher is better)
        'design_to_target_iptm',
        'design_iptm',
        'iptm',
        'protein_iptm',
        # PAE metrics (lower is better)
        'min_interaction_pae',
        'min_design_to_target_pae',
        'interaction_pae',
        # PTM/confidence scores
        'design_ptm',
        'target_ptm',
        'ptm',
        # Designability metrics
        'bb_rmsd',
        'bb_rmsd_design',
        # Liability scores (lower is better)
        'liability_score',
        'liability_num_violations',
        # Additional useful metrics
        'plip_hbonds_refolded',
        'plip_saltbridge_refolded',
        'delta_sasa_refolded',
    ]

    # Keep only columns that exist
    available_columns = [col for col in ranking_columns if col in combined_df.columns]
    ranking_df = combined_df[available_columns].copy()

    return ranking_df


def main():
    print("="*80)
    print("Boltzgen Candidate Scoring (using pre-computed metrics)")
    print("="*80)

    # Load and combine metrics
    df = score_all_directories()

    if df is None or df.empty:
        print("\nNo candidates found to rank.")
        return

    # Sort by key metrics
    # Primary: design_to_target_iptm (higher is better)
    # Secondary: min_interaction_pae (lower is better)
    # Tertiary: liability_score (lower is better)
    df_sorted = df.sort_values(
        by=['design_to_target_iptm', 'min_interaction_pae', 'liability_score'],
        ascending=[False, True, True]
    )

    # Save rankings
    df_sorted.to_csv(OUTPUT_CSV, index=False)

    print(f"\n{'='*80}")
    print(f"SUCCESS! Ranked {len(df)} candidates.")
    print(f"Rankings saved to: {os.path.abspath(OUTPUT_CSV)}")
    print(f"{'='*80}")

    print("\n--- TOP 10 CANDIDATES ---")
    print("\nMetric definitions:")
    print("  design_to_target_iptm: Interface PTM score (0-1, higher = better binding)")
    print("  min_interaction_pae: Min interaction PAE (Å, lower = more confident)")
    print("  liability_score: Developability issues (lower = fewer issues)")
    print()

    # Display top 10 with key metrics
    display_cols = [
        'id',
        'design_to_target_iptm',
        'min_interaction_pae',
        'design_ptm',
        'liability_score',
        'num_design'
    ]
    display_cols = [c for c in display_cols if c in df_sorted.columns]

    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', 30)

    top10 = df_sorted[display_cols].head(10)
    print(top10.to_string(index=False))

    print(f"\n{'='*80}")
    print("\nRanking strategy:")
    print("1. Maximize design_to_target_iptm (binding affinity)")
    print("2. Minimize min_interaction_pae (confidence)")
    print("3. Minimize liability_score (developability)")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    main()
