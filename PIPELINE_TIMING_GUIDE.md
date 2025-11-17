# Pipeline Timing and Logging Guide

## Pipeline Logging

### Log File Location
The pipeline automatically creates a timestamped log file for each run:
```
{output_directory}/pipeline_{YYYYMMDD_HHMMSS}.log
```

Example:
```
results/cxcr4_smart_filter_test/pipeline_20251115_001230.log
```

### Log Contents
The log file contains:
- Timestamp for each major operation
- Stage completion times
- BoltzGen stdout/stderr
- Boltz-2 prediction progress
- Error messages and warnings
- Total pipeline execution time

### Viewing Logs
```bash
# View log file
cat results/your_project/pipeline_*.log

# Follow log in real-time
tail -f results/your_project/pipeline_*.log

# Search for errors
grep "ERROR\|⚠" results/your_project/pipeline_*.log
```

---

## Timing Analysis

### Benchmark Run (H100 GPU)
**Configuration:** 20 designs, 10 candidates, 2 off-targets

| Stage | Time | Notes |
|-------|------|-------|
| BoltzGen Design | 11.6 min | ~35 sec/design |
| BoltzGen Inverse Fold | 1.0 min | ~3 sec/design |
| BoltzGen Folding | 15.7 min | ~47 sec/design |
| BoltzGen Analysis | 0.9 min | - |
| BoltzGen Filtering | 0.2 min | - |
| **Total BoltzGen** | **29.5 min** | **~1.5 min/design** |
| | | |
| Off-Target Predictions | 54 min | ~2.7 min/prediction |
| Selectivity Scoring | 2 min | - |
| **TOTAL PIPELINE** | **~85 min** | **1h 25min** |

---

## Time Estimates for Different Configurations

### Formula
```
Total Time = (num_designs × 1.5) + (num_candidates × num_offtargets × 2.7) + 2 minutes
```

### Common Configurations

| Designs | Candidates | BoltzGen | Off-Target | Total Time |
|---------|-----------|----------|------------|------------|
| 20 | 10 | 30 min | 54 min | **1.4 hours** |
| 50 | 15 | 74 min | 81 min | **2.6 hours** |
| 50 | 20 | 74 min | 108 min | **3.1 hours** |
| 100 | 15 | 148 min | 81 min | **3.8 hours** |
| 100 | 20 | 148 min | 108 min | **4.3 hours** |

*Note: Times assume 2 off-targets and H100 GPU*

---

## Performance Optimization Tips

### 1. Smart Pre-filtering (Default)
Use `num_candidates` to filter designs before off-target scoring:
```yaml
parameters:
  num_designs: 100        # Generate 100 designs
  num_candidates: 15      # Only score top 15 against off-targets
```
**Time saved:** ~70% compared to scoring all designs

### 2. Custom Off-Target Filtering
Fine-tune the number scored against off-targets:
```yaml
parameters:
  num_designs: 100
  num_candidates: 30       # Keep 30 in final output
  num_offtarget_predictions: 15  # Only score top 15 against off-targets
```

### 3. Skip Completed Stages
The pipeline automatically skips:
- BoltzGen if results exist
- Off-target predictions if complete
- Individual predictions if already done

To force re-run, delete the output directory.

---

## Bottleneck Analysis

### Time Distribution
- **Off-target predictions: 60-70%** - Main bottleneck
- **BoltzGen generation: 30-35%** - Secondary bottleneck
- **Scoring & visualization: <5%** - Negligible

### Scaling Behavior

#### Linear Scaling
- `num_designs`: BoltzGen time scales linearly (~1.5 min per design)

#### Multiplicative Scaling (Main Bottleneck!)
- `num_candidates × num_offtargets`: Off-target time scales multiplicatively
  - 10 candidates × 2 targets = 20 predictions × 2.7 min = 54 min
  - 20 candidates × 3 targets = 60 predictions × 2.7 min = 162 min

### Recommendations
1. **Start small:** Test with 20 designs, 10 candidates
2. **Optimize filtering:** Use smart pre-filtering (Option B)
3. **Scale gradually:** Increase to 50-100 designs after validation
4. **Monitor logs:** Watch `pipeline_*.log` for issues

---

## Example Timing Scenarios

### Quick Test (30 minutes)
```yaml
parameters:
  num_designs: 10
  num_candidates: 5
  budget: 10
```
**Time:** ~30 minutes

### Small Production (2.6 hours)
```yaml
parameters:
  num_designs: 50
  num_candidates: 15
  budget: 50
```
**Time:** ~2.6 hours

### Large Production (4.3 hours)
```yaml
parameters:
  num_designs: 100
  num_candidates: 20
  budget: 50
```
**Time:** ~4.3 hours

---

## Monitoring Pipeline Progress

### Real-time Progress
```bash
# Monitor log file
tail -f results/your_project/pipeline_*.log

# Monitor BoltzGen progress
watch -n 5 'ls -lth results/your_project/workbench/intermediate_designs/*.cif | head -20'

# Monitor off-target predictions
watch -n 10 'find results/your_project/workbench/offtarget_predictions -name "confidence_*.json" | wc -l'
```

### Expected Progress Indicators
1. **Stage 3:** Design CIF files appearing in `workbench/intermediate_designs/`
2. **Stage 3.5:** Confidence JSON files in `workbench/offtarget_predictions/*/`
3. **Stage 4:** CSV files in `selectivity_results/`

---

## Troubleshooting

### Pipeline Stuck?
Check the log file:
```bash
tail -100 results/your_project/pipeline_*.log
```

### Out of Memory?
- Reduce `num_designs` or use smaller `budget`
- Check GPU memory: `nvidia-smi`

### Slow Performance?
- Verify GPU is being used (check log for CUDA messages)
- Ensure MSA server is responsive (affects Boltz-2 speed)

---

## Contact & Support
For issues or questions:
- Check log file: `pipeline_*.log`
- Review metrics: `aggregate_metrics_*.csv`
- Inspect failed predictions in `offtarget_predictions/`
