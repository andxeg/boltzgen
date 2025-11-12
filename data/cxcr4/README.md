# CXCR4 Selectivity Test Structures

## Primary Target
- **8u4q_cxcr4.pdb** - CXCR4 receptor (chain R from PDB 8U4Q)
  - Source: REGN7663 Fab bound to CXCR4/Gi complex
  - Use: Primary target for binder generation

## Off-Targets (for selectivity)
- **4mbs_ccr5.pdb** - CCR5 receptor (chain A from PDB 4MBS)
  - Sequence similarity to CXCR4: ~80%
  - Use: Off-target 1 (HIV co-receptor)

- **6lfo_cxcr2.pdb** - CXCR2 receptor (chain R from PDB 6LFO)
  - Sequence similarity to CXCR4: ~60%
  - Use: Off-target 2 (broader negative control)

## Inverse Folding Scaffold
- **8u4q_fab.pdb** - REGN7663 Fab (chains H+L from PDB 8U4Q)
  - Use: Known CXCR4 binder for sequence optimization experiments
  - CDR regions can be redesigned while maintaining framework

## Full Structures (for reference)
- **8u4q_full.pdb** - Complete PDB 8U4Q
- **4mbs_full.pdb** - Complete PDB 4MBS
- **6lfo_full.pdb** - Complete PDB 6LFO
