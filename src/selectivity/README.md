# Selectivity Scorer

Scores binder variants for affinity and selectivity against multiple targets.

## Usage

**Via Pipeline (Recommended):**
```bash
python scripts/optimize_binder.py config.yaml
```

**Standalone:**
```bash
python src/selectivity/scorer.py \
  --variants workbench/designs \
  --primary data/target.pdb \
  --off_targets data/offtarget1.pdb data/offtarget2.pdb \
  --output results/selectivity \
  --num_candidates 30
```

## Output

Creates `aggregate_metrics_selectivity.csv` with:
- `design_id` - Variant identifier
- `sequence` - Amino acid sequence
- `design_to_target_iptm` - Primary target affinity (0-1, higher = better)
- `max_offtarget_iptm` - Worst off-target binding (0-1, lower = better)
- `selectivity_composite` - Overall selectivity score (higher = better)
- `pareto_optimal` - Boolean flag for Pareto frontier membership

## Scoring Formula

```python
# ipTM-based selectivity (higher = more selective)
selectivity_iptm = primary_iptm - max(offtarget_iptms)

# Composite score (weighted average)
selectivity_composite = (
    0.6 * selectivity_iptm +
    0.3 * selectivity_pae +
    0.1 * normalized_sasa
)
```

## How It Works

1. Loads BoltzGen variants from workbench directory
2. Predicts binding to primary target (Boltz-2)
3. Predicts binding to each off-target (Boltz-2)
4. Calculates selectivity scores
5. Ranks designs by composite score
6. Saves top N candidates to CSV

See main README for integration details.
