# Visualization Module

Generates plots and dashboards for selectivity analysis.

## Usage

**Via Pipeline (Recommended):**
```bash
python scripts/optimize_binder.py config.yaml
```
Outputs automatically created in `results/project_name/plots/`

**Standalone:**
```python
from src.visualization.plots import SelectivityVisualizer

visualizer = SelectivityVisualizer(dpi=300)
visualizer.create_comprehensive_dashboard(
    df,
    output_dir="plots",
    title_prefix="Project Name"
)
```

## Generated Plots

**1. Selectivity Dashboard** (`selectivity_dashboard.png`)
- 2x2 comprehensive overview
- Affinity vs selectivity scatter
- Binding profile heatmap
- Distribution and rankings

**2. Pareto Frontier** (`pareto_frontier.png`)
- Multi-objective optimization visualization
- Shows trade-offs between affinity and selectivity

**3. Affinity vs Selectivity** (`affinity_vs_selectivity.png`)
- High-resolution scatter plot
- Color-coded by properties

**4. Binding Profile** (`binding_profile_heatmap.png`)
- Top 30 designs vs all targets
- Red = high binding, dark = low binding

**5. Summary Statistics** (`summary_statistics.csv`)
- Mean, std, min, max, median for key metrics

## Key Features

- Supports both `design_iiptm` and `design_to_target_iptm` column names
- High-resolution output (300 DPI default)
- Publication-ready plots
- Automatic Pareto frontier detection

See main README for complete pipeline documentation.
