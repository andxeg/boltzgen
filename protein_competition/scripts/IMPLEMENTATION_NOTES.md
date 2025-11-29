# Boltz-2 + IPSAE Implementation Notes

## Summary

Successfully created `score_with_boltz_ipsae.py` - a comprehensive pipeline to score Boltzgen candidates using Boltz-2 predictions and IPSAE.

## Key Issues Resolved

### 1. **Boltz Input Format**
- ❌ **Issue**: Boltz doesn't accept CIF files directly
- ✅ **Solution**: Extract sequences from CIF and create YAML input

### 2. **Missing MSAs**
- ❌ **Issue**: `RuntimeError: Missing MSA's in input and --use_msa_server flag not set`
- ✅ **Solution**: Added `--use_msa_server` flag to Boltz command

### 3. **Interactive Prompt in Background**
- ❌ **Issue**: `EOFError` when running script in background
- ✅ **Solution**: Wrapped `input()` in try/except to handle non-interactive mode

## Script Architecture

### Input Processing
```python
CIF file → extract_sequences_from_cif() → {chain_id: sequence}
         → create_boltz_yaml() → YAML file
```

###Boltz Prediction
```python
YAML → boltz predict --use_msa_server
     → predictions/<job_name>/
        ├── <job_name>_model_0.cif
        ├── pae_<job_name>_model_0.npz
        └── confidence_<job_name>_model_0.json
```

### IPSAE Scoring
```python
PAE NPZ + CIF → ipsae.py → parse output → {ipSAE, pDockQ}
```

### Comparison
```python
Boltz scores + Boltzgen CSV → merge → comparison table
```

## Usage Examples

```bash
# Test on 1 candidate
python score_with_boltz_ipsae.py --test-one

# Score top 5 from Boltzgen rankings
python score_with_boltz_ipsae.py --top-n 5

# Score all candidates
python score_with_boltz_ipsae.py --all

# Score specific design
python score_with_boltz_ipsae.py --design-id job_A_medium_00
```

## Performance

- **Per candidate**: ~4-5 minutes
  - MSA generation: 1-2 min
  - Boltz prediction: 2-3 min
  - IPSAE scoring: <10 sec

- **5 candidates**: ~20-25 minutes
- **20 candidates**: ~80-100 minutes

## Output Files

### Main Output
- `boltz_ipsae_scores.csv` - Boltz+IPSAE scores for each candidate
- `boltz_boltzgen_comparison.csv` - Side-by-side comparison

### Comparison Columns
- `ID` - Design identifier
- `Boltz ipSAE` - Interface PTM from Boltz-2 + IPSAE
- `Boltzgen ipTM` - design_to_target_iptm from Boltzgen
- `Boltz pDockQ` - Docking quality from IPSAE
- `Boltzgen PTM` - design_ptm from Boltzgen

## Current Status

✅ Script created with all features
✅ CIF parsing works (extracts chains + sequences correctly)
✅ YAML generation works
✅ Fixed MSA issue (added --use_msa_server)
✅ Fixed interactive prompt issue
🔄 Testing on 1 candidate (in progress)
⏳ Pending: Run on top 5
⏳ Pending: Generate comparison report

## Next Steps

1. **Complete test on 1 candidate** (~5 min)
2. **Verify output quality**
3. **Run on top 5 from Boltzgen rankings**
4. **Analyze correlation between Boltzgen and Boltz+IPSAE metrics**
5. **Generate final comparison report**

## Expected Validation

We expect to see:
- **High correlation** between `design_to_target_iptm` (Boltzgen) and `ipSAE` (Boltz+IPSAE)
- **Similar rankings** for top candidates
- **Slightly different absolute values** (different prediction models)

This will validate that Boltzgen's fast metrics are reliable proxies for the slower Boltz+IPSAE pipeline.
