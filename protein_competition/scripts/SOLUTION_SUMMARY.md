# Scoring Candidates with IPSAE - Solution Summary

## Problem Identified

The original `score_candidates.py` script failed because:

1. **Missing PAE matrices**: Boltzgen's NPZ files (`fold_out_npz/*.npz`) only contain summary statistics like:
   - `interaction_pae` (scalar value)
   - `min_interaction_pae` (scalar value)
   - `design_to_target_iptm` (scalar value)

2. **IPSAE requirements**: The IPSAE script requires full PAE matrices (2D arrays) to compute detailed interface scores, not just summary values.

3. **File format mismatch**: IPSAE expects Boltz-1 format with 'pae' and 'plddt' keys in NPZ files, but Boltzgen uses different key names and doesn't save full matrices.

## Solutions Provided

### ✅ Solution 1: Use Boltzgen Metrics (RECOMMENDED)

**File:** `score_candidates_simple.py`

**What it does:**
- Reads pre-computed metrics from `aggregate_metrics_analyze.csv`
- Ranks candidates by `design_to_target_iptm` (similar to ipSAE)
- Also considers `min_interaction_pae` and `liability_score`
- Instant results - no additional computation

**Usage:**
```bash
cd /lambda/nfs/Nipah-hackathon/andrei/boltzgen/protein_competition/scripts
python score_candidates_simple.py
```

**Output:** `final_rankings_boltzgen.csv`

**Results for your 20 designs:**
- Top candidate: `job_A_medium_00` with design_to_target_iptm = 0.516
- Rankings based on binding affinity, confidence, and developability

**Advantages:**
- ⚡ Fast (instant)
- 📊 Comprehensive metrics already computed
- ✅ No additional dependencies
- 💯 Reliable - uses Boltzgen's native scoring

---

### ⚙️ Solution 2: Full Boltz-2 + IPSAE Pipeline (OPTIONAL)

**File:** `score_candidates.py` (updated)

**What it does:**
1. Takes each CIF file from Boltzgen
2. Runs Boltz-2 prediction to get full PAE matrices
3. Runs IPSAE on the Boltz-2 predictions
4. Generates ipSAE and pDockQ scores

**Changes made:**
- Added `USE_BOLTZ_PREDICTION = True` flag
- Added `run_boltz_prediction()` function
- Updated `score_directory()` to run Boltz-2 before IPSAE
- Added temporary directory management

**Usage:**
```bash
cd /lambda/nfs/Nipah-hackathon/andrei/boltzgen/protein_competition/scripts
python score_candidates.py
```

**⚠️ Warning:** This is VERY SLOW (~3-5 minutes per structure × 20 = 60-100 minutes total)

**When to use:**
- You specifically need ipSAE/pDockQ scores for publication
- You want to compare with other methods that use IPSAE
- You need detailed interface analysis with PyMOL scripts

---

## Comparison of Metrics

| What you want to measure | Boltzgen Metric | IPSAE Metric |
|--------------------------|----------------|--------------|
| Binding affinity | `design_to_target_iptm` | `ipSAE` |
| Interface confidence | `min_interaction_pae` | PAE cutoff analysis |
| Complex quality | `design_ptm`, `iptm` | `pDockQ` |
| Designability | `bb_rmsd`, `design_ptm` | N/A |
| Developability | `liability_score` | N/A |

**Both are measuring essentially the same thing** - the interface PTM (ipTM) and PAE scores.

---

## Recommended Workflow

### For quick screening (RECOMMENDED):
```bash
python score_candidates_simple.py
```

### For detailed IPSAE analysis (if needed):
```bash
# Test on one candidate first
python test_score_one.py

# If successful, run on all
python score_candidates.py
```

---

## Files Created/Modified

1. ✅ `score_candidates_simple.py` - NEW: Fast scoring using Boltzgen metrics
2. ✅ `score_candidates.py` - UPDATED: Added Boltz-2 prediction support
3. ✅ `test_score_one.py` - NEW: Test Boltz+IPSAE on single candidate
4. ✅ `README_SCORING.md` - NEW: Documentation
5. ✅ `SOLUTION_SUMMARY.md` - NEW: This file

---

## Next Steps

1. **Already done**: `score_candidates_simple.py` has been run
   - Output: `final_rankings_boltzgen.csv`
   - Top candidate: job_A_medium_00

2. **Optional**: Run `test_score_one.py` to verify Boltz+IPSAE works
   - Currently running in background

3. **If needed**: Run full IPSAE scoring with `score_candidates.py`
   - Only if you specifically need ipSAE scores
   - Will take 1-2 hours for all 20 candidates

---

## Quick Check: See Your Top Candidates

```bash
# View top 10 ranked candidates
head -11 final_rankings_boltzgen.csv | column -t -s','

# Or get just the IDs and scores
cut -d',' -f1,5 final_rankings_boltzgen.csv | head -11
```
