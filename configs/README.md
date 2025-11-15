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
  num_candidates: 30       # Top N for final output
  num_offtarget_predictions: 30  # OPTIONAL: How many to score against off-targets
                           # - If not set: uses num_candidates (default, recommended)
                           # - If "all": scores all num_designs (exhaustive)
                           # - If number: scores that many top candidates
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

## Off-Target Prediction Strategies

The `num_offtarget_predictions` parameter controls how many designs are scored against off-targets (Stage 3.5). This is separate from `num_candidates` which controls the final output size.

**Option A: Score All Designs (Exhaustive)**
```yaml
parameters:
  num_designs: 20
  num_candidates: 5
  num_offtarget_predictions: "all"  # Score all 20 designs
```
- Use when: You want selectivity data for every design
- Performance: Slowest (all designs × off-targets predictions)
- Example: `examples/option_a_score_all.yaml`

**Option B: Smart Pre-filtering (Default, Recommended)**
```yaml
parameters:
  num_designs: 100
  num_candidates: 10
  # num_offtarget_predictions: NOT SET (defaults to num_candidates)
```
- Use when: You want efficient computation (default)
- Performance: Fast (only top candidates scored)
- Rationale: Poor primary binders won't be selected anyway
- Example: `examples/option_b_smart_filter.yaml`

**Option C: Custom Filtering**
```yaml
parameters:
  num_designs: 50
  num_candidates: 15
  num_offtarget_predictions: 20  # Score top 20
```
- Use when: You want fine-tuned control
- Performance: Configurable trade-off
- Example: `examples/option_custom.yaml`

## Available Configs

**Production Configs:**
- `cxcr4_existing_binder_full_test.yaml` - CXCR4 optimization (recommended)
- `cxcr4_denovo_example.yaml` - De novo design example
- `regn7663_optimization_small.yaml` - Small test run

**Strategy Examples:**
- `examples/option_a_score_all.yaml` - Score all designs (exhaustive)
- `examples/option_b_smart_filter.yaml` - Smart pre-filtering (default)
- `examples/option_custom.yaml` - Custom configuration

See main README for detailed documentation.
