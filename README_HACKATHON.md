# 🎯 Fine-Grained Binder Design Platform

**Production-ready pipeline for therapeutic binder optimization with selectivity guidance**

---

## 📋 Table of Contents

1. [Quick Start](#quick-start)
2. [Pipeline Overview](#pipeline-overview)
3. [Installation & Setup](#installation--setup)
4. [Running the Pipelines](#running-the-pipelines)
5. [Configuration Guide](#configuration-guide)
6. [Understanding Outputs](#understanding-outputs)
7. [Parameter Tuning](#parameter-tuning)
8. [Troubleshooting](#troubleshooting)
9. [Examples](#examples)

---

## 🚀 Quick Start

### Run Existing Binder Optimization (RECOMMENDED)

```bash
# Single command execution
./optimize_binder.sh configs/cxcr4_optimization_example.yaml
```

### Run De Novo Binder Design

```bash
./optimize_binder.sh configs/cxcr4_denovo_example.yaml
```

**Expected Time**: 20-40 minutes for 100 designs

---

## 🔬 Pipeline Overview

### Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    YAML Configuration File                       │
│  (targets, scaffold, parameters, scoring weights, output)       │
└──────────────────────┬──────────────────────────────────────────┘
                       │
                       ▼
         ┌─────────────────────────────┐
         │   [1] Input Validation      │
         │   • Check PDB files exist   │
         │   • Validate parameters     │
         │   • Verify design regions   │
         └──────────┬──────────────────┘
                    │
                    ▼
         ┌─────────────────────────────┐
         │   [2] BoltzGen Generation   │
         │   ┌─────────────────────┐   │
         │   │ Existing Binder     │   │  (Inverse Folding on Scaffold)
         │   │  - Scaffold PDB     │   │
         │   │  - CDR regions      │   │
         │   │  - Target PDB       │   │
         │   └─────────────────────┘   │
         │           OR                │
         │   ┌─────────────────────┐   │
         │   │   De Novo Design    │   │  (Generate from Scratch)
         │   │  - Binder length    │   │
         │   │  - Target PDB       │   │
         │   └─────────────────────┘   │
         │                             │
         │   Output: 100+ designs      │
         └──────────┬──────────────────┘
                    │
                    ▼
         ┌─────────────────────────────┐
         │  [3] Selectivity Scoring    │
         │   • Primary target binding  │
         │   • Off-target scoring      │
         │   • Composite metrics       │
         │   • Select top N            │
         └──────────┬──────────────────┘
                    │
                    ▼
         ┌─────────────────────────────┐
         │  [4] Pareto Optimization    │
         │   • Find Pareto frontier    │
         │   • Multi-layer ranking     │
         │   • Trade-off analysis      │
         │   • Diversity selection     │
         └──────────┬──────────────────┘
                    │
                    ▼
         ┌─────────────────────────────┐
         │  [5] Output Generation      │
         │   • ranked_designs.csv      │
         │   • results.json            │
         │   • plots/                  │
         │     - selectivity_dashboard │
         │     - pareto_frontier       │
         │     - binding_heatmap       │
         └─────────────────────────────┘
```

### Two Modes

**Mode 1: Existing Binder Optimization** (Priority)
```
Input: Known binder scaffold + target
Process: Inverse folding on CDR regions
Output: Optimized variants with improved affinity/selectivity
Use case: Improve existing antibody, nanobody, or DARPin
```

**Mode 2: De Novo Binder Design** (Extension)
```
Input: Target only
Process: Generate new binders from scratch
Output: Novel binders ranked by affinity/selectivity
Use case: Design completely new therapeutic binder
```

---

## 🛠 Installation & Setup

### Prerequisites

1. **BoltzGen** (already installed in your environment)
2. **Python 3.8+** with packages:
   - pandas
   - numpy
   - matplotlib
   - seaborn
   - pyyaml

### Verify Installation

```bash
# Check BoltzGen
boltzgen --version

# Check Python packages
python -c "import pandas, numpy, matplotlib, seaborn, yaml; print('✓ All packages installed')"
```

### Directory Structure

```
boltzgen/
├── optimize_binder.sh           # Main entry point (bash wrapper)
├── scripts/
│   └── optimize_binder.py       # Pipeline orchestrator
├── src/
│   ├── selectivity/
│   │   └── scorer.py            # Selectivity scoring
│   ├── optimization/
│   │   └── pareto.py            # Pareto optimization
│   ├── output/
│   │   └── json_generator.py   # JSON output generation
│   └── visualization/
│       └── plots.py             # Visualization dashboard
├── configs/
│   ├── cxcr4_optimization_example.yaml  # Existing binder example
│   └── cxcr4_denovo_example.yaml        # De novo example
├── data/
│   └── cxcr4/
│       ├── 8u4q_cxcr4.pdb      # Primary target
│       ├── 8u4q_fab.pdb        # Scaffold
│       ├── 4mbs_ccr5.pdb       # Off-target 1
│       └── 6lfo_cxcr2.pdb      # Off-target 2
└── results/                     # Output directory (created automatically)
```

---

## 🏃 Running the Pipelines

### Method 1: Using the Bash Wrapper (RECOMMENDED)

```bash
./optimize_binder.sh configs/your_config.yaml
```

### Method 2: Direct Python Execution

```bash
python scripts/optimize_binder.py configs/your_config.yaml
```

### Monitoring Progress

The pipeline provides real-time feedback:

```
================================================================================
Fine-Grained Binder Design Pipeline
================================================================================
Project: CXCR4_optimization
Type: existing_binder
Output: results/cxcr4_optimization
================================================================================

[1/5] Validating input files...
  ✓ All input files found

[2/5] Generating BoltzGen configuration...
  ✓ BoltzGen config saved: results/cxcr4_optimization/boltzgen_config.yaml

[3/5] Running BoltzGen design generation...
  Generating 100 designs with budget 50...
  This may take 20-40 minutes...
  ✓ BoltzGen generation complete

[4/5] Running selectivity scoring...
  Scoring top 30 candidates...
  Primary target: CXCR4
  Off-targets: CCR5, CXCR2
  ✓ Selectivity scoring complete

[5/5] Generating outputs...
  ✓ CSV results: results/cxcr4_optimization/ranked_designs.csv
  ✓ Plots: results/cxcr4_optimization/plots/
  ✓ All outputs saved to: results/cxcr4_optimization

================================================================================
✓ PIPELINE COMPLETE!
================================================================================
```

---

## ⚙️ Configuration Guide

### YAML Configuration Structure

```yaml
# PROJECT METADATA
project:
  name: "CXCR4_optimization"       # Project identifier
  type: "existing_binder"          # or "de_novo"

# TARGET SPECIFICATION
targets:
  # Positive targets (want to bind strongly)
  positive:
    - name: "CXCR4"                # Display name
      pdb: "data/cxcr4/8u4q_cxcr4.pdb"
      chain: "R"

  # Negative targets (want to avoid binding)
  negative:
    - name: "CCR5"
      pdb: "data/cxcr4/4mbs_ccr5.pdb"
      chain: "A"
    - name: "CXCR2"
      pdb: "data/cxcr4/6lfo_cxcr2.pdb"
      chain: "A"

# SCAFFOLD CONFIGURATION (for existing_binder mode only)
scaffold:
  pdb: "data/cxcr4/8u4q_fab.pdb"
  chains: ["H", "L"]              # Heavy and light chains
  design_regions:                 # CDR regions to optimize
    - chain: "H"
      residues: "26-35,50-65,95-102"  # CDR H1,H2,H3
    - chain: "L"
      residues: "24-34,50-56,89-97"   # CDR L1,L2,L3

# BINDER SPECIFICATION (for de_novo mode only)
binder_spec:
  length: "110-130"               # Length range for nanobody
  type: "nanobody"                # Binder type

# DESIGN PARAMETERS
parameters:
  num_designs: 100                # Total designs to generate (50-2000)
  num_candidates: 30              # Top N for selectivity scoring (10-100)
  budget: 50                      # BoltzGen sampling budget (20-100)

# SCORING CONFIGURATION
scoring:
  weights:
    affinity: 0.6                 # Primary target binding weight (0.0-1.0)
    selectivity: 0.3              # Off-target avoidance weight (0.0-1.0)
    properties: 0.1               # Biophysical properties weight (0.0-1.0)
  multi_objective: true           # Enable Pareto ranking

# OUTPUT CONFIGURATION
output:
  directory: "results/cxcr4_optimization"
  format: ["csv", "json", "plots"]
```

### Required vs Optional Fields

**Required for ALL modes**:
- `project.name`
- `project.type`
- `targets.positive` (at least 1)
- `parameters.num_designs`
- `output.directory`

**Required for `existing_binder` mode**:
- `scaffold.pdb`
- `scaffold.chains`
- `scaffold.design_regions`

**Required for `de_novo` mode**:
- `binder_spec.length`
- `binder_spec.type`

**Optional**:
- `targets.negative` (selectivity scoring disabled if not provided)
- `parameters.num_candidates` (default: 30)
- `parameters.budget` (default: 50)
- `scoring.*` (uses defaults if not provided)

---

## 📊 Understanding Outputs

### Output Directory Structure

```
results/cxcr4_optimization/
├── ranked_designs.csv              # Main results table
├── results.json                    # Structured JSON for web UI
├── plots/
│   ├── selectivity_dashboard.png   # 4-panel overview
│   ├── pareto_frontier_primary.png # Affinity vs selectivity trade-offs
│   ├── affinity_vs_selectivity.png # Scatter plot
│   └── binding_profile_heatmap.png # Multi-target binding matrix
└── workbench/                      # Intermediate files (BoltzGen output)
```

### CSV Output (`ranked_designs.csv`)

Key columns:
- **design_id**: Unique identifier
- **sequence**: Full amino acid sequence
- **design_to_target_iptm**: Primary target affinity (0-1, higher = better)
- **min_design_to_target_pae**: Primary target confidence (Å, lower = better)
- **max_offtarget_iptm**: Worst off-target binding (0-1, lower = better)
- **selectivity_composite**: Overall selectivity score (higher = better)
- **pareto_optimal**: Boolean flag for Pareto frontier membership
- **pareto_rank**: Layer number (0 = frontier, 1 = second layer, etc.)

### JSON Output (`results.json`)

Structure:
```json
{
  "project": {
    "name": "CXCR4_optimization",
    "timestamp": "2025-11-11T20:00:00Z",
    "targets": {...}
  },
  "summary": {
    "total_designs": 100,
    "pareto_optimal": 12,
    "mean_affinity": 0.65,
    "mean_selectivity": 0.42
  },
  "ranked_designs": [
    {
      "rank": 1,
      "design_id": "design_003",
      "sequence": "DVV...",
      "metrics": {
        "affinity": {"primary_iptm": 0.75, "primary_pae": 5.2},
        "selectivity": {"CCR5_iptm": 0.32, "selectivity_score": 0.48},
        "properties": {"sasa": 3850}
      },
      "composite_score": 0.68,
      "pareto_optimal": true
    },
    ...
  ],
  "pareto_analysis": {
    "frontier_designs": [1, 3, 7, 12],
    "trade_offs": {...}
  }
}
```

### Visualization Plots

**selectivity_dashboard.png** (4-panel overview):
- Top-left: Affinity vs selectivity scatter
- Top-right: Binding profile heatmap
- Bottom-left: Selectivity distribution
- Bottom-right: Ranking strategy comparison

**pareto_frontier_primary.png**:
- Shows trade-off between affinity and selectivity
- Red points = Pareto-optimal designs
- Red dashed line = Pareto frontier

---

## 🎛 Parameter Tuning

### Balancing Speed vs Quality

| Parameter | Fast | Balanced | High Quality |
|-----------|------|----------|--------------|
| num_designs | 50 | 100-200 | 500-2000 |
| budget | 20 | 50 | 100 |
| num_candidates | 10 | 30 | 100 |
| **Time** | ~10 min | ~30 min | ~2-6 hours |

### Adjusting Scoring Weights

**High Affinity Priority** (strong binding to target):
```yaml
scoring:
  weights:
    affinity: 0.8      # Prioritize binding strength
    selectivity: 0.1   # Less concern for off-targets
    properties: 0.1
```

**High Selectivity Priority** (avoid off-targets):
```yaml
scoring:
  weights:
    affinity: 0.4
    selectivity: 0.5   # Prioritize selectivity
    properties: 0.1
```

**Balanced** (default):
```yaml
scoring:
  weights:
    affinity: 0.6
    selectivity: 0.3
    properties: 0.1
```

### Design Region Selection (Existing Binder Mode)

**Conservative** (small changes):
```yaml
design_regions:
  - chain: "H"
    residues: "95-102"  # Only CDR H3 (most variable)
```

**Moderate** (typical):
```yaml
design_regions:
  - chain: "H"
    residues: "26-35,50-65,95-102"  # All heavy chain CDRs
  - chain: "L"
    residues: "24-34,50-56,89-97"   # All light chain CDRs
```

**Aggressive** (extensive redesign):
```yaml
design_regions:
  - chain: "H"
    residues: "20-110"  # Entire variable domain
  - chain: "L"
    residues: "20-110"
```

---

## 🔧 Troubleshooting

### Common Issues

**Issue 1: "Missing input files"**
```
ERROR: Missing input files:
  - data/cxcr4/8u4q_cxcr4.pdb
```
**Solution**: Verify all PDB paths in your YAML config are correct and files exist.

**Issue 2: "BoltzGen command not found"**
```
boltzgen: command not found
```
**Solution**: Activate your conda/virtual environment with BoltzGen installed.

**Issue 3: "No designs generated"**
```
ERROR: No designs found in workbench directory
```
**Solution**: Check BoltzGen logs. May need to adjust parameters or check GPU availability.

**Issue 4: "Zero affinity scores"**
```
Mean affinity: 0.000
```
**Solution**: This is an old issue - should be fixed. If you still see this, verify you're using the latest version of `src/selectivity/scorer.py` with `design_to_target_iptm` metric.

### Debug Mode

Enable verbose logging:
```bash
python scripts/optimize_binder.py configs/your_config.yaml --verbose
```

Check intermediate files:
```bash
ls -la results/your_project/workbench/
```

---

## 📝 Examples

### Example 1: Quick Test (10 minutes)

```yaml
# config/quick_test.yaml
project:
  name: "quick_test"
  type: "de_novo"

targets:
  positive:
    - name: "CXCR4"
      pdb: "data/cxcr4/8u4q_cxcr4.pdb"
      chain: "R"

binder_spec:
  length: "110-130"
  type: "nanobody"

parameters:
  num_designs: 20        # Small for quick testing
  num_candidates: 5
  budget: 20

output:
  directory: "results/quick_test"
  format: ["csv", "plots"]
```

Run:
```bash
./optimize_binder.sh configs/quick_test.yaml
```

### Example 2: High-Quality Production Run

```yaml
# config/production.yaml
project:
  name: "CXCR4_production"
  type: "existing_binder"

targets:
  positive:
    - name: "CXCR4"
      pdb: "data/cxcr4/8u4q_cxcr4.pdb"
      chain: "R"
  negative:
    - name: "CCR5"
      pdb: "data/cxcr4/4mbs_ccr5.pdb"
      chain: "A"
    - name: "CXCR2"
      pdb: "data/cxcr4/6lfo_cxcr2.pdb"
      chain: "A"

scaffold:
  pdb: "data/cxcr4/8u4q_fab.pdb"
  chains: ["H", "L"]
  design_regions:
    - chain: "H"
      residues: "26-35,50-65,95-102"
    - chain: "L"
      residues: "24-34,50-56,89-97"

parameters:
  num_designs: 500       # High quality
  num_candidates: 100
  budget: 100

scoring:
  weights:
    affinity: 0.6
    selectivity: 0.3
    properties: 0.1
  multi_objective: true

output:
  directory: "results/cxcr4_production"
  format: ["csv", "json", "plots"]
```

### Example 3: Swapping Targets (Easy Portability)

To optimize for a different target, just change the PDB files:

```yaml
targets:
  positive:
    - name: "HER2"
      pdb: "data/her2/6og2_her2.pdb"  # Change this
      chain: "A"
  negative:
    - name: "EGFR"
      pdb: "data/egfr/1nql_egfr.pdb"  # Change this
      chain: "A"
```

Everything else stays the same!

---

## 🎯 Best Practices

### For Quick Testing

1. **Start with quick test** (20 designs) to verify pipeline
2. **Run production** (100-200 designs) for real results
3. **Focus on visualizations** - dashboard provides comprehensive overview
4. **Highlight Pareto frontier** - shows multi-objective optimization
5. **Demonstrate portability** - swap targets easily

### For Production Use

1. **Validate inputs** carefully (correct chains, residues)
2. **Start conservative** on design regions
3. **Iterate** on scoring weights based on results
4. **Use selectivity** for any therapeutic application
5. **Analyze Pareto frontier** - don't just pick rank 1

### Performance Tips

- Use **GPU** for BoltzGen (much faster)
- **Parallel execution** if multiple targets (future feature)
- **Cache results** - reuse BoltzGen output for different scoring weights

---

## 📚 Additional Resources

- **BoltzGen Documentation**: https://github.com/jwohlwend/boltz
- **Pipeline Architecture**: See diagram in this README (Section: Pipeline Overview)

---

## 🆘 Getting Help

**Issues? Questions?**
1. Check this README first
2. Review implementation summary
3. Check error messages carefully
4. Examine intermediate outputs in `workbench/`

**Support**:
- All components tested and working
- Example configs provided
- Full pipeline validated

---

## ✅ Pre-Flight Checklist

Before running the pipeline:

- [ ] BoltzGen installed and accessible
- [ ] Python packages installed (pandas, numpy, matplotlib, seaborn, yaml)
- [ ] PDB files available in correct locations
- [ ] YAML config created with correct paths
- [ ] Output directory specified and writable
- [ ] GPU available (for faster execution)

---

**Ready to go! Run your first pipeline:**

```bash
./optimize_binder.sh configs/cxcr4_optimization_example.yaml
```

**Happy optimizing! 🚀**
