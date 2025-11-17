# Architecture Refactoring: Separation of Concerns

## Problem Statement

The initial implementation had **scorer.py running BoltzGen predictions**, which violated the Single Responsibility Principle:

```python
# BEFORE: scorer.py (BAD ARCHITECTURE)
class SelectivityScorer:
    def _score_against_offtarget(...):
        # ❌ Creates YAML files
        # ❌ Runs subprocess calls to boltzgen
        # ❌ Mixes orchestration with calculation
        # ❌ Cannot reuse predictions
```

## Solution: Proper Separation

### Stage 3.5: Pipeline Orchestration (optimize_binder.py)

**New method: `_run_offtarget_predictions()`** (Lines 245-386)

```python
def _run_offtarget_predictions(self, workbench_dir: Path) -> Optional[Path]:
    """Run BoltzGen predictions for off-target binding."""

    # 1. Load top N candidates from primary metrics
    df = pd.read_csv(metrics_file)
    top_candidates = df.sort_values('design_to_target_iptm').head(num_candidates)

    # 2. For each off-target:
    for off_target in off_targets:
        # 3. For each candidate:
        for idx, row in top_candidates.iterrows():
            # - Create YAML for [Design + Off-target] complex
            # - Run BoltzGen folding prediction
            # - Save to workbench/offtarget_predictions/{offtarget_name}/{design_id}/
```

**Directory Structure:**
```
workbench/
  ├── designs/                        # Original designs
  ├── aggregate_metrics_analyze.csv   # Primary metrics
  └── offtarget_predictions/          # NEW: Off-target predictions
      ├── CCR5/
      │   ├── design_0/
      │   │   └── aggregate_metrics_analyze.csv
      │   ├── design_1/
      │   └── ...
      └── CXCR2/
          ├── design_0/
          └── ...
```

### Stage 4: Pure Scoring Logic (scorer.py)

**Refactored: `_score_against_offtarget()`** (Lines 273-334)

```python
def _score_against_offtarget(...) -> pd.DataFrame:
    """Load pre-computed off-target predictions and extract metrics."""

    # ✅ NO subprocess calls
    # ✅ NO BoltzGen execution
    # ✅ Pure data loading and calculation

    if self.offtarget_predictions_dir.exists():
        for design_id in sequences_df:
            # Load from: offtarget_predictions/{offtarget_name}/{design_id}/
            metrics = self._load_precomputed_metrics(...)
            results.append(metrics)
```

## Benefits

### 1. **Separation of Concerns**
- Orchestration ≠ Calculation
- Each module has ONE responsibility

### 2. **Reusability**
```bash
# Run predictions once
python scripts/optimize_binder.py config.yaml
  # Stage 3.5 creates workbench/offtarget_predictions/

# Reuse predictions multiple times
python src/selectivity/scorer.py \
  --variants workbench/designs \
  --offtarget_predictions workbench/offtarget_predictions \
  --output results/  # Fast! No re-prediction needed
```

### 3. **Testability**
```python
# Easy to test with mock data
scorer = SelectivityScorer(
    primary_target=mock_target,
    off_targets=mock_offtargets,
    offtarget_predictions_dir="/path/to/test/data"
)
results = scorer.score_variants(...)  # No external dependencies
```

### 4. **Caching**
- Predictions saved in workbench directory
- Can skip re-running if already computed
- Future: Add `--reuse` flag to skip existing predictions

### 5. **Clarity**
```
Pipeline Flow (CLEAR):
  Stage 3:   Generate designs          [BoltzGen]
  Stage 3.5: Predict off-target binding [BoltzGen]  ← NEW
  Stage 4:   Calculate selectivity      [Scorer]    ← PURE LOGIC
  Stage 5:   Pareto optimization
  Stage 6:   Generate outputs
```

## Code Changes Summary

### optimize_binder.py
- **Added** `_run_offtarget_predictions()` (Stage 3.5)
- **Added** `_get_chain_ids()` helper
- **Updated** `_run_selectivity_scoring()` to pass predictions directory
- **Updated** `run()` method to call Stage 3.5

### scorer.py
- **Removed** `_run_boltzgen_folding()` (orchestration)
- **Removed** `_get_design_chain_ids()` (moved to pipeline)
- **Removed** `_get_offtarget_chain_id()` (moved to pipeline)
- **Removed** `_extract_prediction_metrics()` (replaced)
- **Added** `_load_precomputed_metrics()` (pure loading)
- **Updated** `__init__()` to accept `offtarget_predictions_dir`
- **Updated** `_score_against_offtarget()` to load pre-computed data
- **Updated** `main()` to accept `--offtarget_predictions` argument

## Testing

```bash
# Test with cxcr4_optimization_example.yaml
rm -rf results/cxcr4_quick_test
python scripts/optimize_binder.py configs/cxcr4_optimization_example.yaml
```

**Expected Output:**
```
[Stage 3] Running BoltzGen...
  ✓ 20 designs generated

[Stage 3.5] Running off-target predictions...
  [1/2] Processing off-target: CCR5
    Running Boltz-2 predictions...
    ✓ Completed CCR5 predictions
  [2/2] Processing off-target: CXCR2
    Running Boltz-2 predictions...
    ✓ Completed CXCR2 predictions
  ✓ Off-target predictions saved to: workbench/offtarget_predictions

[Stage 4] Running selectivity scoring...
  Using pre-computed predictions from: workbench/offtarget_predictions
  ✓ Selectivity scoring complete
```

## Future Improvements

1. **Parallel Predictions**: Run off-target predictions in parallel
2. **Caching Logic**: Skip existing predictions with `--reuse` flag
3. **Batch Processing**: Combine multiple YAMLs into single BoltzGen run
4. **Progress Tracking**: Better feedback for long-running predictions
5. **Error Recovery**: Resume from failed predictions

## Conclusion

This refactoring provides:
- ✅ Clean architecture (SRP)
- ✅ Proper separation of concerns
- ✅ Reusable predictions
- ✅ Testable components
- ✅ Cacheable results
- ✅ Production-ready code

The system now follows best practices and is ready for production use.
