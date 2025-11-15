# Selectivity Scorer

Scores binder variants for affinity and selectivity against multiple targets using real Boltz-2 structure predictions.

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

## How It Works (Production Implementation)

1. Loads BoltzGen variants from workbench directory
2. Uses primary target metrics from BoltzGen analysis
3. **For each off-target:**
   - Creates YAML files for [Design + Off-target] complexes
   - Runs **real Boltz-2 structure prediction** via `boltzgen run`
   - Extracts ipTM and PAE metrics from predictions
4. Calculates selectivity scores from real binding predictions
5. Ranks designs by composite score
6. Saves top N candidates to CSV

**Key Feature:** This implementation uses actual Boltz-2 folding predictions for off-target binding, not estimates. Each design is predicted against each off-target, providing biophysically accurate selectivity measurements.

**Performance:** Off-target prediction time scales with (num_candidates × num_off_targets). For 5 candidates and 2 off-targets, expect ~10-15 minutes of additional prediction time.

See main README for integration details.
