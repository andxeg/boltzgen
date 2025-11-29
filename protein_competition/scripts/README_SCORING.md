# Candidate Scoring Methods

This directory contains two approaches for scoring Boltzgen candidates:

## Method 1: Use Pre-computed Boltzgen Metrics (RECOMMENDED)

**Script:** `score_candidates_simple.py`

**Advantages:**
- ✅ Fast - no additional computation needed
- ✅ Uses Boltzgen's native metrics computed during generation
- ✅ Includes comprehensive metrics: affinity, confidence, designability, liability
- ✅ Already includes `design_to_target_iptm` which is similar to ipSAE

**Usage:**
```bash
python score_candidates_simple.py
```

**Key Metrics:**
- `design_to_target_iptm`: Interface PTM score (0-1, higher = better binding affinity)
- `min_interaction_pae`: Minimum interaction PAE in Angstroms (lower = more confident)
- `design_ptm`: Design confidence (0-1, higher = better)
- `liability_score`: Developability issues (lower = fewer problems)

**Output:** `final_rankings_boltzgen.csv`

---

## Method 2: Use Boltz-2 + IPSAE (SLOW, for detailed analysis)

**Script:** `score_candidates.py`

**Advantages:**
- Provides ipSAE and pDockQ scores specifically
- Can be used for detailed interface analysis
- Generates PyMOL scripts for visualization

**Disadvantages:**
- ⚠️ Very slow - runs Boltz-2 prediction on each candidate (~2-5 min per structure)
- ⚠️ Requires GPU
- ⚠️ May take hours for 20+ candidates

**Usage:**
```bash
# Edit score_candidates.py and set USE_BOLTZ_PREDICTION = True
python score_candidates.py
```

This will:
1. Take each CIF file from Boltzgen
2. Run Boltz-2 prediction to get full PAE matrices
3. Run IPSAE on the predictions
4. Generate rankings with ipSAE and pDockQ scores

**Output:** `final_rankings.csv`

---

## Comparison: Boltzgen Metrics vs IPSAE

| Metric | Boltzgen | IPSAE |
|--------|----------|-------|
| Binding affinity | `design_to_target_iptm` | `ipSAE` |
| Interface confidence | `min_interaction_pae` | PAE-based |
| Docking quality | `design_ptm` | `pDockQ` |
| Computation time | Instant (pre-computed) | Hours (requires Boltz-2) |

**Bottom line:** For most purposes, the Boltzgen metrics are sufficient and much faster. Use IPSAE only if you need the specific scores for comparison or publication.

---

## Current Issue with score_candidates.py

The original `score_candidates.py` script had an issue:
- Boltzgen NPZ files don't contain full PAE matrices (only summary statistics)
- IPSAE requires full PAE matrices
- Solution: Run Boltz-2 predictions first (set `USE_BOLTZ_PREDICTION = True`)

This has been fixed in the updated version.
