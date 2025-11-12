"""
Pareto Optimization for Multi-Objective Binder Design
======================================================

Implements Pareto frontier analysis for balancing multiple objectives:
- Affinity (primary target binding)
- Selectivity (off-target avoidance)
- Biophysical properties (solubility, stability, etc.)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional
import matplotlib.pyplot as plt


class ParetoOptimizer:
    """
    Multi-objective Pareto optimization for binder design.

    Identifies non-dominated solutions and provides trade-off analysis.
    """

    def __init__(
        self,
        objectives: Optional[List[str]] = None,
        maximize: Optional[Dict[str, bool]] = None
    ):
        """
        Initialize Pareto optimizer.

        Args:
            objectives: List of objective column names
                       Default: ['design_to_target_iptm', 'selectivity_composite']
            maximize: Dict mapping objective names to True (maximize) or False (minimize)
        """
        self.objectives = objectives or [
            'design_to_target_iptm',      # Higher is better
            'selectivity_composite'       # Higher is better
        ]

        self.maximize = maximize or {
            'design_to_target_iptm': True,
            'selectivity_composite': True,
            'min_design_to_target_pae': False,  # Lower is better
            'delta_sasa_refolded': True,        # Higher buried surface
            'solubility': True,
            'stability': True
        }

    def find_pareto_frontier(
        self,
        df: pd.DataFrame,
        objectives: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Identify Pareto-optimal (non-dominated) designs.

        A design is Pareto-optimal if no other design is better in all objectives.

        Args:
            df: DataFrame with design metrics
            objectives: Objective columns (uses self.objectives if None)

        Returns:
            DataFrame with 'pareto_optimal' column added
        """
        if objectives is None:
            objectives = [obj for obj in self.objectives if obj in df.columns]

        if len(objectives) < 2:
            raise ValueError("Need at least 2 objectives for Pareto optimization")

        print(f"\nPareto Optimization:")
        print(f"  Objectives: {', '.join(objectives)}")

        # Prepare objective matrix
        obj_matrix = df[objectives].values

        # Flip sign for minimization objectives
        for i, obj in enumerate(objectives):
            if not self.maximize.get(obj, True):
                obj_matrix[:, i] = -obj_matrix[:, i]

        # Find Pareto frontier
        pareto_mask = self._is_pareto_efficient(obj_matrix)

        df['pareto_optimal'] = pareto_mask
        num_pareto = pareto_mask.sum()

        print(f"  ✓ Found {num_pareto} Pareto-optimal designs ({num_pareto/len(df)*100:.1f}%)")

        return df

    def _is_pareto_efficient(self, costs: np.ndarray) -> np.ndarray:
        """
        Find Pareto-efficient points (maximization).

        Args:
            costs: Objective matrix (n_designs x n_objectives)
                  All objectives should be maximization (flip sign if minimizing)

        Returns:
            Boolean mask of Pareto-efficient points
        """
        is_efficient = np.ones(costs.shape[0], dtype=bool)

        for i, c in enumerate(costs):
            if is_efficient[i]:
                # Mark dominated points as not efficient
                # A point is dominated if another point is >= in all objectives
                # and strictly > in at least one
                is_efficient[is_efficient] = np.any(
                    costs[is_efficient] > c, axis=1
                ) | np.all(costs[is_efficient] == c, axis=1)

        return is_efficient

    def rank_pareto_layers(
        self,
        df: pd.DataFrame,
        objectives: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Rank designs by Pareto layers (frontier depth).

        Layer 0: Pareto frontier
        Layer 1: Pareto frontier after removing layer 0
        etc.

        Args:
            df: DataFrame with design metrics
            objectives: Objective columns

        Returns:
            DataFrame with 'pareto_rank' column
        """
        if objectives is None:
            objectives = [obj for obj in self.objectives if obj in df.columns]

        df_copy = df.copy()
        df_copy['pareto_rank'] = -1

        layer = 0
        remaining_mask = np.ones(len(df_copy), dtype=bool)

        while remaining_mask.sum() > 0:
            # Find Pareto frontier of remaining points
            remaining_df = df_copy[remaining_mask]

            obj_matrix = remaining_df[objectives].values
            for i, obj in enumerate(objectives):
                if not self.maximize.get(obj, True):
                    obj_matrix[:, i] = -obj_matrix[:, i]

            pareto_mask = self._is_pareto_efficient(obj_matrix)

            # Assign rank
            pareto_indices = remaining_df.index[pareto_mask]
            df_copy.loc[pareto_indices, 'pareto_rank'] = layer

            # Remove current frontier
            remaining_mask[pareto_indices] = False
            layer += 1

            if layer > 10:  # Safety limit
                break

        print(f"  ✓ Ranked into {layer} Pareto layers")

        return df_copy

    def analyze_trade_offs(
        self,
        df: pd.DataFrame,
        objectives: Optional[List[str]] = None
    ) -> Dict:
        """
        Analyze trade-offs between objectives on Pareto frontier.

        Args:
            df: DataFrame with 'pareto_optimal' column
            objectives: Objective columns

        Returns:
            Dict with trade-off analysis
        """
        if objectives is None:
            objectives = [obj for obj in self.objectives if obj in df.columns]

        pareto_df = df[df['pareto_optimal']]

        analysis = {
            'num_pareto': len(pareto_df),
            'num_total': len(df),
            'pareto_fraction': len(pareto_df) / len(df),
            'objective_ranges': {},
            'correlations': {},
            'trade_offs': {}
        }

        # Objective ranges on Pareto frontier
        for obj in objectives:
            analysis['objective_ranges'][obj] = {
                'min': float(pareto_df[obj].min()),
                'max': float(pareto_df[obj].max()),
                'mean': float(pareto_df[obj].mean()),
                'std': float(pareto_df[obj].std())
            }

        # Pairwise correlations on Pareto frontier
        for i, obj1 in enumerate(objectives):
            for obj2 in objectives[i+1:]:
                corr = pareto_df[[obj1, obj2]].corr().iloc[0, 1]
                analysis['correlations'][f'{obj1}_vs_{obj2}'] = float(corr)

                # Classify trade-off
                if corr < -0.3:
                    trade_off = "Strong trade-off (negative correlation)"
                elif corr > 0.3:
                    trade_off = "Aligned objectives (positive correlation)"
                else:
                    trade_off = "Independent objectives"

                analysis['trade_offs'][f'{obj1}_vs_{obj2}'] = trade_off

        return analysis

    def plot_pareto_frontier_2d(
        self,
        df: pd.DataFrame,
        obj_x: str,
        obj_y: str,
        output_path: Optional[Path] = None
    ):
        """
        Visualize 2D Pareto frontier.

        Args:
            df: DataFrame with 'pareto_optimal' column
            obj_x: X-axis objective
            obj_y: Y-axis objective
            output_path: Save path (shows plot if None)
        """
        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot all designs
        non_pareto = df[~df['pareto_optimal']]
        ax.scatter(
            non_pareto[obj_x],
            non_pareto[obj_y],
            c='lightgray',
            s=50,
            alpha=0.6,
            label='Non-Pareto designs'
        )

        # Plot Pareto frontier
        pareto = df[df['pareto_optimal']]
        ax.scatter(
            pareto[obj_x],
            pareto[obj_y],
            c='red',
            s=100,
            alpha=0.8,
            edgecolors='darkred',
            linewidth=2,
            label=f'Pareto optimal ({len(pareto)} designs)',
            zorder=10
        )

        # Connect Pareto points
        if len(pareto) > 1:
            pareto_sorted = pareto.sort_values(obj_x)
            ax.plot(
                pareto_sorted[obj_x],
                pareto_sorted[obj_y],
                'r--',
                alpha=0.5,
                linewidth=1,
                label='Pareto frontier',
                zorder=5
            )

        ax.set_xlabel(obj_x, fontsize=12)
        ax.set_ylabel(obj_y, fontsize=12)
        ax.set_title(f'Pareto Frontier: {obj_x} vs {obj_y}', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"  ✓ Saved Pareto plot: {output_path}")
            plt.close()
        else:
            plt.show()

    def select_representative_designs(
        self,
        df: pd.DataFrame,
        num_designs: int = 5,
        objectives: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Select representative designs from Pareto frontier.

        Uses k-means clustering in objective space to find diverse solutions.

        Args:
            df: DataFrame with 'pareto_optimal' column
            num_designs: Number of representative designs to select
            objectives: Objective columns

        Returns:
            DataFrame with top representative designs
        """
        if objectives is None:
            objectives = [obj for obj in self.objectives if obj in df.columns]

        pareto_df = df[df['pareto_optimal']].copy()

        if len(pareto_df) <= num_designs:
            return pareto_df

        # Normalize objectives to [0, 1]
        obj_matrix = pareto_df[objectives].values
        obj_normalized = (obj_matrix - obj_matrix.min(axis=0)) / (
            obj_matrix.max(axis=0) - obj_matrix.min(axis=0) + 1e-10
        )

        # Simple diversity selection: maximize distance in objective space
        selected_indices = []

        # Start with best overall design (highest sum of normalized objectives)
        scores = obj_normalized.sum(axis=1)
        selected_indices.append(scores.argmax())

        # Iteratively select design furthest from existing selections
        for _ in range(num_designs - 1):
            if len(selected_indices) >= len(pareto_df):
                break

            selected_points = obj_normalized[selected_indices]

            # Find design with maximum minimum distance to selected points
            distances = np.min(
                np.linalg.norm(
                    obj_normalized[:, None, :] - selected_points[None, :, :],
                    axis=2
                ),
                axis=1
            )

            # Don't select already selected designs
            distances[selected_indices] = -1

            next_idx = distances.argmax()
            selected_indices.append(next_idx)

        return pareto_df.iloc[selected_indices]

    def optimize(
        self,
        df: pd.DataFrame,
        objectives: Optional[List[str]] = None,
        output_dir: Optional[Path] = None
    ) -> Dict:
        """
        Complete Pareto optimization workflow.

        Args:
            df: DataFrame with design metrics
            objectives: Objective columns
            output_dir: Directory to save results and plots

        Returns:
            Dict with optimization results
        """
        print("\n" + "="*60)
        print("Multi-Objective Pareto Optimization")
        print("="*60)

        if objectives is None:
            objectives = [obj for obj in self.objectives if obj in df.columns]

        # Find Pareto frontier
        df = self.find_pareto_frontier(df, objectives)

        # Rank by Pareto layers
        df = self.rank_pareto_layers(df, objectives)

        # Analyze trade-offs
        analysis = self.analyze_trade_offs(df, objectives)

        print(f"\nTrade-off Analysis:")
        for trade_off, desc in analysis['trade_offs'].items():
            corr = analysis['correlations'][trade_off]
            print(f"  {trade_off}: {desc} (r={corr:.2f})")

        # Select representative designs
        representatives = self.select_representative_designs(df, num_designs=5, objectives=objectives)

        print(f"\nRepresentative Designs:")
        for i, (idx, row) in enumerate(representatives.iterrows(), 1):
            obj_str = ', '.join([f"{obj}={row[obj]:.3f}" for obj in objectives[:3]])
            print(f"  {i}. {row['design_id']}: {obj_str}")

        # Generate visualizations
        if output_dir and len(objectives) >= 2:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            # Plot primary trade-off (first 2 objectives)
            self.plot_pareto_frontier_2d(
                df,
                objectives[0],
                objectives[1],
                output_dir / 'pareto_frontier_primary.png'
            )

        results = {
            'df': df,
            'analysis': analysis,
            'representatives': representatives,
            'num_pareto': len(df[df['pareto_optimal']]),
            'pareto_indices': df[df['pareto_optimal']].index.tolist()
        }

        print("="*60 + "\n")

        return results


if __name__ == '__main__':
    # Example usage
    import argparse

    parser = argparse.ArgumentParser(description='Pareto optimization for binder designs')
    parser.add_argument('--input', type=Path, required=True, help='Input CSV with metrics')
    parser.add_argument('--output', type=Path, required=True, help='Output directory')
    parser.add_argument(
        '--objectives',
        nargs='+',
        default=['design_to_target_iptm', 'selectivity_composite'],
        help='Objective columns'
    )

    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.input)

    # Run optimization
    optimizer = ParetoOptimizer(objectives=args.objectives)
    results = optimizer.optimize(df, output_dir=args.output)

    # Save results
    args.output.mkdir(parents=True, exist_ok=True)
    results['df'].to_csv(args.output / 'pareto_optimized.csv', index=False)
    results['representatives'].to_csv(args.output / 'pareto_representatives.csv', index=False)

    print(f"✓ Results saved to {args.output}")
