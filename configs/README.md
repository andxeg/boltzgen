# Configuration Files

YAML configurations for binder design pipelines.

## Quick Start

**Existing Binder Optimization:**
```bash
python scripts/optimize_binder.py configs/cxcr4_existing_binder_full_test.yaml
```

**De Novo Design:**
```bash
python scripts/optimize_binder.py configs/cxcr4_denovo_example.yaml
```

## Configuration Structure

```yaml
project:
  name: "project_name"
  type: "existing_binder"  # or "de_novo"

targets:
  positive:                # Targets to bind
    - name: "CXCR4"
      pdb: "data/cxcr4/8u4q_cxcr4.pdb"
  negative:                # Off-targets to avoid
    - name: "CCR5"
      pdb: "data/cxcr4/4mbs_ccr5.pdb"

scaffold:                  # For existing_binder only
  pdb: "data/cxcr4/8u4q_fab.pdb"
  design_regions:
    - chain: "H"
      residues: "26-35,50-65,95-102"

parameters:
  num_designs: 100         # Total designs to generate
  num_candidates: 30       # Top N for selectivity scoring
  budget: 50               # BoltzGen sampling budget

scoring:
  multi_objective: true    # Enable Pareto optimization
  weights:
    affinity: 0.6
    selectivity: 0.3
    properties: 0.1

output:
  directory: "results/project_name"
  format: ["csv", "json", "plots"]
```

## Available Configs

- `cxcr4_existing_binder_full_test.yaml` - CXCR4 optimization (recommended)
- `cxcr4_denovo_example.yaml` - De novo design example
- `regn7663_optimization_small.yaml` - Small test run

See main README for detailed documentation.
