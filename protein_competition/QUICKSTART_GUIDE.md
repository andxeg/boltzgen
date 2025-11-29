# Boltzgen → Scoring Pipeline - Quick Start Guide

## Step-by-Step Workflow

### Step 1: Generate Candidates with Boltzgen

```bash
# Navigate to your project directory
cd /path/to/your/project

# Create your design specification YAML
# See example: boltzgen/example/vanilla_peptide_with_target_binding_site/beetletert.yaml
nano my_design.yaml

# Run Boltzgen to generate candidates
boltzgen predict my_design.yaml \
  --output results/my_experiment \
  --num-designs 20 \
  --inverse-folding-budget 5

# This creates: results/my_experiment/
#   ├── intermediate_designs/          # Initial designs
#   ├── intermediate_designs_inverse_folded/  # Refined designs (use these!)
#   │   ├── refold_cif/                # Structure files
#   │   ├── fold_out_npz/              # Confidence metrics
#   │   └── aggregate_metrics_analyze.csv  # All metrics
#   └── final_ranked_designs/          # Top candidates
```

---

### Step 2: Score Candidates

**Choose ONE of the following methods:**

---

## Method A: FAST Scoring (Recommended for screening)

**Uses:** Pre-computed Boltzgen metrics
**Time:** < 1 second
**Best for:** Quick screening, initial ranking

```bash
# Navigate to scripts directory
cd results/my_experiment
mkdir -p scripts
cd scripts

# Copy the scoring script
cp /path/to/score_candidates_simple.py .

# Edit configuration
nano score_candidates_simple.py
# Set: RESULT_DIRS = ["../"]

# Run scoring
python score_candidates_simple.py

# View results
cat final_rankings_boltzgen.csv
```

**Output:**
- `final_rankings_boltzgen.csv` - All candidates ranked by:
  - `design_to_target_iptm` - Binding affinity (higher = better)
  - `min_interaction_pae` - Confidence (lower = better)
  - `liability_score` - Developability (lower = better)

---

## Method B: DETAILED Scoring (For validation/publication)

**Uses:** Boltz-2 predictions + IPSAE scoring
**Time:** ~4-5 minutes per candidate
**Best for:** Final validation, detailed analysis, publication

```bash
# Navigate to scripts directory
cd results/my_experiment/scripts

# Copy the Boltz+IPSAE script
cp /path/to/score_with_boltz_ipsae.py .

# Option 1: Test on 1 candidate first
python score_with_boltz_ipsae.py --test-one

# Option 2: Score top N from Boltzgen rankings
python score_with_boltz_ipsae.py --top-n 5

# Option 3: Score all candidates (takes hours!)
python score_with_boltz_ipsae.py --all

# Option 4: Score specific design
python score_with_boltz_ipsae.py --design-id my_design_00

# View results
cat boltz_ipsae_scores.csv
cat boltz_boltzgen_comparison.csv
```

**Output:**
- `boltz_ipsae_scores.csv` - Boltz+IPSAE scores
- `boltz_boltzgen_comparison.csv` - Side-by-side comparison

---

## Method C: BOTH (Best Practice)

**Recommended workflow:**

```bash
# 1. Fast screening with Boltzgen metrics
python score_candidates_simple.py

# 2. Validate top 5-10 with Boltz+IPSAE
python score_with_boltz_ipsae.py --top-n 10

# 3. Compare results
cat boltz_boltzgen_comparison.csv
```

---

## Complete Example

```bash
# Full pipeline from scratch
cd /lambda/nfs/Nipah-hackathon/andrei/boltzgen

# 1. GENERATE (example already exists)
# boltzgen predict example/vanilla_peptide/design.yaml --output results/test

# 2. FAST SCORING
cd protein_competition/scripts

# Edit RESULT_DIRS in script
echo 'RESULT_DIRS = ["../results/20designs_5budget"]' > config.txt

python score_candidates_simple.py
# ✅ Output: final_rankings_boltzgen.csv (instant)

# 3. DETAILED SCORING (top 5)
python score_with_boltz_ipsae.py --top-n 5
# ⏱️ Takes ~20-25 minutes
# ✅ Output: boltz_ipsae_scores.csv, boltz_boltzgen_comparison.csv

# 4. VIEW RESULTS
echo "=== Boltzgen Rankings ==="
head -6 final_rankings_boltzgen.csv | column -t -s','

echo "=== Boltz+IPSAE Comparison ==="
cat boltz_boltzgen_comparison.csv | column -t -s','
```

---

## Quick Reference

### Scoring Metrics Explained

| Metric | Source | Range | Better | Meaning |
|--------|--------|-------|--------|---------|
| `design_to_target_iptm` | Boltzgen | 0-1 | Higher | Binding affinity |
| `min_interaction_pae` | Boltzgen | Å | Lower | Prediction confidence |
| `design_ptm` | Boltzgen | 0-1 | Higher | Structure quality |
| `liability_score` | Boltzgen | 0+ | Lower | Developability issues |
| `ipSAE` | Boltz+IPSAE | 0-1 | Higher | Interface quality |
| `pDockQ` | Boltz+IPSAE | 0-1 | Higher | Docking quality |

### File Locations

```
your_project/
├── design.yaml              # Input specification
├── results/
│   └── experiment_name/
│       ├── intermediate_designs_inverse_folded/
│       │   ├── refold_cif/  # ← Structure files for scoring
│       │   └── aggregate_metrics_analyze.csv  # ← Boltzgen metrics
│       └── scripts/
│           ├── score_candidates_simple.py
│           ├── score_with_boltz_ipsae.py
│           ├── final_rankings_boltzgen.csv      # ← Fast results
│           └── boltz_ipsae_scores.csv           # ← Detailed results
```

---

## Troubleshooting

### "No such file or directory: IPSAE"
```bash
# Check IPSAE installation
ls /lambda/nfs/Nipah-hackathon/andrei/IPSAE/ipsae.py

# Update IPSAE_SCRIPT path in score_with_boltz_ipsae.py
```

### "Boltz prediction failed"
```bash
# Ensure boltz is installed and accessible
which boltz
conda list | grep boltz

# Check GPU availability
nvidia-smi
```

### "Missing rankings CSV"
```bash
# Run fast scoring first
python score_candidates_simple.py

# This creates final_rankings_boltzgen.csv needed for --top-n
```

---

## Performance Tips

1. **Use fast scoring first** - Get instant results, identify top candidates
2. **Validate selectively** - Only run Boltz+IPSAE on top 5-10
3. **Run in parallel** - Score multiple experiments simultaneously
4. **Keep temp files** - Useful for debugging (answer 'y' when prompted)

---

## Script Parameters

### score_candidates_simple.py
```python
# Edit these in the script
RESULT_DIRS = ["path/to/results"]  # Your experiment directories
OUTPUT_CSV = "final_rankings_boltzgen.csv"  # Output filename
```

### score_with_boltz_ipsae.py
```python
# Command line options
--test-one           # Test on 1 candidate
--top-n N            # Score top N from Boltzgen rankings
--all                # Score all candidates
--design-id ID       # Score specific design

# Edit in script if needed
PAE_CUTOFF = 10
DIST_CUTOFF = 10
BOLTZ_DIFFUSION_SAMPLES = 1  # Increase for better quality (slower)
```

---

## Expected Runtime

| Task | Time | Notes |
|------|------|-------|
| Boltzgen generation (20 designs) | 30-60 min | Depends on design size |
| Fast scoring (any number) | < 1 sec | Uses pre-computed metrics |
| Boltz+IPSAE (1 candidate) | 4-5 min | MSA + prediction + IPSAE |
| Boltz+IPSAE (top 5) | 20-25 min | Can run in background |
| Boltz+IPSAE (20 candidates) | 80-100 min | Run overnight |

---

**That's it!** You now have a complete pipeline from generation to scoring. Start with fast scoring, then validate top candidates with Boltz+IPSAE if needed.
