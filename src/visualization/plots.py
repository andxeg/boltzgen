#!/usr/bin/env python3
"""
Visualization Module for Selective Binder Design

Creates comprehensive visualizations for multi-target selectivity analysis:
1. Affinity vs Selectivity scatter plots
2. Binding profile heatmaps
3. Selectivity distribution
4. Pareto frontier visualization
5. Property comparison charts
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List
import argparse


class SelectivityVisualizer:
    """
    Create publication-quality visualizations for selectivity analysis.

    Visualizes:
    - Trade-offs between affinity and selectivity
    - Binding profiles across targets
    - Multi-objective optimization results
    - Comparison to baseline (original binder)
    """

    def __init__(self, figsize: tuple = (16, 12), dpi: int = 300):
        """
        Initialize visualizer.

        Args:
            figsize: Default figure size (width, height)
            dpi: Resolution for saved figures
        """
        self.figsize = figsize
        self.dpi = dpi

        # Set publication-quality style
        plt.style.use('seaborn-v0_8-darkgrid')
        sns.set_palette("husl")

    def create_comprehensive_dashboard(
        self,
        metrics_df: pd.DataFrame,
        output_dir: Path,
        baseline_metrics: Optional[pd.DataFrame] = None,
        title_prefix: str = "REGN7663 Optimization"
    ):
        """
        Create comprehensive 2x2 dashboard with all key visualizations.

        Args:
            metrics_df: DataFrame with selectivity metrics
            output_dir: Directory to save plots
            baseline_metrics: Optional baseline (original binder) metrics
            title_prefix: Prefix for plot titles
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)

        print("\n" + "="*60)
        print("Creating Selectivity Visualizations")
        print("="*60)

        # Create 2x2 subplot layout
        fig, axes = plt.subplots(2, 2, figsize=self.figsize)
        fig.suptitle(f'{title_prefix}: Multi-Target Selectivity Analysis',
                    fontsize=16, fontweight='bold', y=0.995)

        # 1. Top-left: Affinity vs Selectivity scatter
        print("[1/4] Affinity vs Selectivity scatter...")
        self._plot_affinity_vs_selectivity(metrics_df, axes[0, 0], baseline_metrics)

        # 2. Top-right: Binding profile heatmap
        print("[2/4] Binding profile heatmap...")
        self._plot_binding_heatmap(metrics_df, axes[0, 1])

        # 3. Bottom-left: Selectivity distribution
        print("[3/4] Selectivity distribution...")
        self._plot_selectivity_distribution(metrics_df, axes[1, 0])

        # 4. Bottom-right: Property comparison
        print("[4/4] Property comparison...")
        self._plot_property_comparison(metrics_df, axes[1, 1])

        plt.tight_layout()

        # Save dashboard
        dashboard_file = output_dir / "selectivity_dashboard.png"
        plt.savefig(dashboard_file, dpi=self.dpi, bbox_inches='tight')
        print(f"\n✓ Dashboard saved: {dashboard_file}")

        plt.close()

        # Create individual high-res plots
        print("\nCreating individual plots...")
        self._create_individual_plots(metrics_df, output_dir, baseline_metrics)

        print("\n" + "="*60)
        print("Visualization complete!")
        print("="*60 + "\n")

    def _plot_affinity_vs_selectivity(
        self,
        df: pd.DataFrame,
        ax: plt.Axes,
        baseline: Optional[pd.DataFrame] = None
    ):
        """Plot affinity vs selectivity with color-coded properties."""

        # Check for affinity column (support both naming conventions)
        affinity_col = None
        if 'design_to_target_iptm' in df.columns:
            affinity_col = 'design_to_target_iptm'
        elif 'design_iiptm' in df.columns:
            affinity_col = 'design_iiptm'

        # Check required columns
        if affinity_col is None or 'selectivity_composite' not in df.columns:
            ax.text(0.5, 0.5, 'Missing data', ha='center', va='center', transform=ax.transAxes)
            return

        # Color by solubility if available
        if 'solubility' in df.columns:
            scatter = ax.scatter(
                df[affinity_col],
                df['selectivity_composite'],
                c=df['solubility'],
                cmap='viridis',
                alpha=0.7,
                s=50,
                edgecolors='black',
                linewidth=0.5
            )
            plt.colorbar(scatter, ax=ax, label='Solubility')
        else:
            ax.scatter(
                df[affinity_col],
                df['selectivity_composite'],
                alpha=0.7,
                s=50,
                edgecolors='black',
                linewidth=0.5
            )

        # Add baseline if provided
        if baseline is not None:
            baseline_affinity_col = None
            if 'design_to_target_iptm' in baseline.columns:
                baseline_affinity_col = 'design_to_target_iptm'
            elif 'design_iiptm' in baseline.columns:
                baseline_affinity_col = 'design_iiptm'

            if baseline_affinity_col:
                ax.scatter(
                    baseline[baseline_affinity_col].mean(),
                    0,  # Baseline has no selectivity optimization
                    color='red',
                    s=200,
                    marker='*',
                    edgecolors='black',
                    linewidth=2,
                    label='Original REGN7663',
                    zorder=10
                )

        # Reference lines
        ax.axhline(0, color='red', linestyle='--', linewidth=2, alpha=0.5,
                   label='No selectivity')
        ax.axhline(df['selectivity_composite'].median(), color='green',
                   linestyle='--', linewidth=1, alpha=0.5, label='Median selectivity')

        ax.set_xlabel('Primary Target Affinity (ipTM)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Selectivity Score', fontsize=11, fontweight='bold')
        ax.set_title('Affinity vs Selectivity', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)

    def _plot_binding_heatmap(self, df: pd.DataFrame, ax: plt.Axes):
        """Plot binding profile heatmap for top designs."""

        # Get top 30 by selectivity
        top_df = df.nlargest(30, 'selectivity_composite')

        # Find target columns (support both naming conventions)
        affinity_col = 'design_to_target_iptm' if 'design_to_target_iptm' in df.columns else 'design_iiptm'
        target_cols = [affinity_col] if affinity_col in df.columns else []
        offtarget_cols = [c for c in df.columns if 'offtarget_iptm' in c]

        all_cols = target_cols + offtarget_cols

        if len(all_cols) < 2:
            ax.text(0.5, 0.5, 'Insufficient target data', ha='center', va='center',
                   transform=ax.transAxes)
            return

        # Prepare heatmap data
        heatmap_data = top_df[all_cols].T

        # Rename rows for clarity
        row_labels = ['CXCR4 (Primary)'] + [col.replace('offtarget_iptm_', '').replace('_', ' ')
                                             for col in offtarget_cols]

        # Plot heatmap
        sns.heatmap(
            heatmap_data,
            ax=ax,
            cmap='RdYlGn',
            cbar_kws={'label': 'ipTM'},
            yticklabels=row_labels,
            xticklabels=False,
            vmin=0.3,
            vmax=0.9,
            linewidths=0.5,
            linecolor='gray'
        )

        ax.set_title('Top 30 Designs: Binding Profile', fontsize=12, fontweight='bold')
        ax.set_xlabel('Design Variants (ranked by selectivity)', fontsize=10)
        ax.set_ylabel('Target', fontsize=11, fontweight='bold')

    def _plot_selectivity_distribution(self, df: pd.DataFrame, ax: plt.Axes):
        """Plot selectivity score distribution."""

        if 'selectivity_composite' not in df.columns:
            ax.text(0.5, 0.5, 'Missing selectivity data', ha='center', va='center',
                   transform=ax.transAxes)
            return

        # Histogram
        ax.hist(
            df['selectivity_composite'],
            bins=30,
            alpha=0.7,
            color='steelblue',
            edgecolor='black',
            linewidth=1.2
        )

        # Reference lines
        ax.axvline(0, color='red', linestyle='--', linewidth=2, label='No selectivity')
        ax.axvline(
            df['selectivity_composite'].mean(),
            color='orange',
            linestyle='-',
            linewidth=2,
            label=f'Mean: {df["selectivity_composite"].mean():.3f}'
        )
        ax.axvline(
            df['selectivity_composite'].median(),
            color='green',
            linestyle='-',
            linewidth=2,
            label=f'Median: {df["selectivity_composite"].median():.3f}'
        )

        ax.set_xlabel('Selectivity Score', fontsize=11, fontweight='bold')
        ax.set_ylabel('Number of Designs', fontsize=11, fontweight='bold')
        ax.set_title('Selectivity Distribution', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=9)
        ax.grid(axis='y', alpha=0.3)

    def _plot_property_comparison(self, df: pd.DataFrame, ax: plt.Axes):
        """Compare properties: selective ranking vs affinity ranking."""

        # Find affinity column (support both naming conventions)
        affinity_col = 'design_to_target_iptm' if 'design_to_target_iptm' in df.columns else 'design_iiptm'

        if affinity_col not in df.columns:
            ax.text(0.5, 0.5, 'Missing affinity data', ha='center', va='center',
                   transform=ax.transAxes)
            return

        # Top 30 by selectivity
        top_selective = df.nlargest(30, 'selectivity_composite')

        # Top 30 by affinity (no selectivity consideration)
        top_affinity = df.nlargest(30, affinity_col)

        # Comparison metrics
        metrics = {}

        if affinity_col in df.columns:
            metrics['Primary\nAffinity'] = [
                top_selective[affinity_col].mean(),
                top_affinity[affinity_col].mean()
            ]

        if 'max_offtarget_iptm' in df.columns:
            metrics['Max\nOff-Target'] = [
                top_selective['max_offtarget_iptm'].mean(),
                top_affinity['max_offtarget_iptm'].mean()
            ]

        if 'selectivity_composite' in df.columns:
            metrics['Selectivity\nScore'] = [
                top_selective['selectivity_composite'].mean(),
                top_affinity['selectivity_composite'].mean()
            ]

        if 'solubility' in df.columns:
            metrics['Solubility'] = [
                top_selective['solubility'].mean(),
                top_affinity['solubility'].mean()
            ]

        if not metrics:
            ax.text(0.5, 0.5, 'Insufficient metric data', ha='center', va='center',
                   transform=ax.transAxes)
            return

        # Create grouped bar chart
        comparison_df = pd.DataFrame(
            metrics,
            index=['Selective Ranking', 'Affinity Ranking']
        )

        comparison_df.T.plot(
            kind='bar',
            ax=ax,
            color=['steelblue', 'coral'],
            alpha=0.8,
            edgecolor='black',
            linewidth=1.2
        )

        ax.set_ylabel('Score', fontsize=11, fontweight='bold')
        ax.set_title('Ranking Strategy Comparison', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=9)
        ax.grid(axis='y', alpha=0.3)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    def _create_individual_plots(
        self,
        df: pd.DataFrame,
        output_dir: Path,
        baseline: Optional[pd.DataFrame] = None
    ):
        """Create individual high-resolution plots."""

        # 1. Large affinity vs selectivity plot
        fig, ax = plt.subplots(figsize=(10, 8))
        self._plot_affinity_vs_selectivity(df, ax, baseline)
        plt.tight_layout()
        plt.savefig(output_dir / "affinity_vs_selectivity.png", dpi=self.dpi, bbox_inches='tight')
        plt.close()
        print("  ✓ affinity_vs_selectivity.png")

        # 2. Large heatmap
        fig, ax = plt.subplots(figsize=(12, 6))
        self._plot_binding_heatmap(df, ax)
        plt.tight_layout()
        plt.savefig(output_dir / "binding_profile_heatmap.png", dpi=self.dpi, bbox_inches='tight')
        plt.close()
        print("  ✓ binding_profile_heatmap.png")

        # 3. Pareto frontier (if multiple objectives)
        affinity_col = 'design_to_target_iptm' if 'design_to_target_iptm' in df.columns else 'design_iiptm'
        if affinity_col in df.columns and 'selectivity_composite' in df.columns:
            fig, ax = plt.subplots(figsize=(10, 8))
            self._plot_pareto_frontier(df, ax)
            plt.tight_layout()
            plt.savefig(output_dir / "pareto_frontier.png", dpi=self.dpi, bbox_inches='tight')
            plt.close()
            print("  ✓ pareto_frontier.png")

    def _plot_pareto_frontier(self, df: pd.DataFrame, ax: plt.Axes):
        """Plot Pareto frontier for multi-objective optimization."""

        # Find affinity column (support both naming conventions)
        affinity_col = 'design_to_target_iptm' if 'design_to_target_iptm' in df.columns else 'design_iiptm'

        # For 2D: affinity vs selectivity
        x = df[affinity_col].values
        y = df['selectivity_composite'].values

        # Find Pareto frontier
        pareto_mask = np.ones(len(df), dtype=bool)
        for i in range(len(df)):
            for j in range(len(df)):
                if i != j:
                    # j dominates i if j is better in both objectives
                    if x[j] >= x[i] and y[j] >= y[i] and (x[j] > x[i] or y[j] > y[i]):
                        pareto_mask[i] = False
                        break

        # Plot all points
        ax.scatter(x, y, alpha=0.4, s=30, label='All designs', color='gray')

        # Highlight Pareto optimal
        ax.scatter(
            x[pareto_mask],
            y[pareto_mask],
            s=100,
            color='red',
            edgecolors='black',
            linewidth=2,
            label=f'Pareto optimal ({pareto_mask.sum()} designs)',
            zorder=10
        )

        # Connect Pareto frontier
        pareto_x = x[pareto_mask]
        pareto_y = y[pareto_mask]
        sorted_indices = np.argsort(pareto_x)
        ax.plot(
            pareto_x[sorted_indices],
            pareto_y[sorted_indices],
            'r--',
            linewidth=2,
            alpha=0.6,
            label='Pareto frontier'
        )

        ax.set_xlabel('Primary Affinity (ipTM)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Selectivity Score', fontsize=12, fontweight='bold')
        ax.set_title('Pareto Frontier: Affinity vs Selectivity', fontsize=14, fontweight='bold')
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)

    def create_summary_table(self, df: pd.DataFrame, output_dir: Path):
        """Create summary statistics table."""

        output_dir = Path(output_dir)

        summary = {
            'Metric': [],
            'Mean': [],
            'Std': [],
            'Min': [],
            'Max': [],
            'Median': []
        }

        # Support both column naming conventions
        affinity_col = 'design_to_target_iptm' if 'design_to_target_iptm' in df.columns else 'design_iiptm'

        metrics_to_summarize = [
            (affinity_col, 'Primary Affinity (ipTM)'),
            ('max_offtarget_iptm', 'Max Off-Target Affinity'),
            ('selectivity_composite', 'Selectivity Score'),
            ('solubility', 'Solubility'),
            ('mean_plddt_design', 'Mean pLDDT')
        ]

        for col, label in metrics_to_summarize:
            if col in df.columns:
                summary['Metric'].append(label)
                summary['Mean'].append(f"{df[col].mean():.3f}")
                summary['Std'].append(f"{df[col].std():.3f}")
                summary['Min'].append(f"{df[col].min():.3f}")
                summary['Max'].append(f"{df[col].max():.3f}")
                summary['Median'].append(f"{df[col].median():.3f}")

        summary_df = pd.DataFrame(summary)

        # Save to CSV
        summary_file = output_dir / "summary_statistics.csv"
        summary_df.to_csv(summary_file, index=False)

        # Print to console
        print("\n" + "="*60)
        print("Summary Statistics")
        print("="*60)
        print(summary_df.to_string(index=False))
        print("="*60 + "\n")

        return summary_df


def main():
    """Command-line interface for visualization."""

    parser = argparse.ArgumentParser(
        description="Create visualizations for selectivity analysis"
    )
    parser.add_argument(
        '--input', type=Path, required=True,
        help="Selectivity metrics CSV file (aggregate_metrics_selectivity.csv)"
    )
    parser.add_argument(
        '--output', type=Path, required=True,
        help="Output directory for plots"
    )
    parser.add_argument(
        '--baseline', type=Path,
        help="Optional baseline metrics CSV for comparison"
    )
    parser.add_argument(
        '--title', type=str, default="REGN7663 Optimization",
        help="Title prefix for plots"
    )
    parser.add_argument(
        '--dpi', type=int, default=300,
        help="Resolution for saved figures (default: 300)"
    )

    args = parser.parse_args()

    # Load data
    print(f"Loading data from: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  ✓ Loaded {len(df)} variants\n")

    baseline = None
    if args.baseline:
        print(f"Loading baseline from: {args.baseline}")
        baseline = pd.read_csv(args.baseline)
        print(f"  ✓ Loaded baseline\n")

    # Create visualizer
    visualizer = SelectivityVisualizer(dpi=args.dpi)

    # Generate visualizations
    visualizer.create_comprehensive_dashboard(
        df,
        args.output,
        baseline,
        args.title
    )

    # Create summary table
    visualizer.create_summary_table(df, args.output)

    print(f"\n✓ All visualizations saved to: {args.output}")


if __name__ == '__main__':
    main()
