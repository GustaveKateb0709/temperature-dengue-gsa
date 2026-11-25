# Temperature-Forced Dengue Transmission Model with LHS-PRCC Global Sensitivity Analysis

Reproducible code accompanying the manuscript *Temperature-driven dengue
transmission: a reproducible modelling and global sensitivity analysis
framework with a ranking-matrix protocol* (in preparation).

The repository contains:

- A temperature-forced compartmental dengue transmission model (human and
  mosquito hosts, seasonal temperature forcing, annual importation pulse).
- A Latin Hypercube Sampling - Partial Rank Correlation Coefficient (LHS-PRCC)
  global sensitivity analysis pipeline around R0 and dynamic epidemic outputs.
- A **ranking-matrix protocol study**: a systematic comparison of how PRCC
  parameter rankings change under alternative input-range designs (published
  ranges vs equal-width ranges vs +/-20 % ranges), across static-R0 and
  dynamic-output responses.
- A verification test suite that pins every headline number reported in the
  manuscript.

All results are generated with fixed seeds; re-running the pipeline
reproduces the manuscript numbers exactly.

## Repository layout

```
.
├── src/                    # Core model and main analysis pipeline
│   ├── config.py           # Global configuration (seed, paths, defaults)
│   ├── model.py            # Temperature-forced dengue compartmental model
│   ├── r0.py               # R0 computation and temperature response curve
│   ├── sensitivity.py      # LHS-PRCC sensitivity analysis
│   ├── scenarios.py        # Epidemic scenario definitions
│   ├── run_all.py          # Main pipeline: runs all analyses, 5 figures
│   ├── make_figures.py     # Figure generation (300 dpi PNG + PDF)
│   ├── figure_export.py    # Self-contained Matplotlib export helper
│   └── style_presets.py    # Shared figure style presets
├── tests/
│   └── test_verification.py  # Assertion suite pinning manuscript numbers
├── new_analysis/           # Dynamic-output sensitivity and ranking-matrix study
│   ├── dv_common.py        # Shared helpers (model wrapper, LHS, JSON I/O)
│   ├── make_static_reference.py   # -> static_prcc_reference.json
│   ├── run_dynamic_prcc.py        # -> dynamic_prcc.json
│   ├── run_scaling.py             # -> nv_scaling.json, nv_scaling.md
│   ├── run_clamp.py               # -> clamp_sensitivity.json
│   ├── run_equalwidth.py          # -> equalwidth_prcc.json
│   ├── run_equalwidth50.py        # -> equalwidth50_prcc.json
│   ├── run_static_equalwidth.py   # -> static_equalwidth_prcc.json
│   ├── run_constant.py            # -> constant_param_prcc.json
│   ├── run_rank_matrix.py         # -> rank_matrix.json, rank_matrix_summary.md
│   ├── make_summary.py            # -> SUMMARY.md
│   ├── stamp_src_hashes.py        # Adds source hashes / provenance to results
│   ├── *.json                     # Archived result files (as used in the paper)
│   ├── SUMMARY.md                 # Human-readable summary of all analyses
│   └── rank_matrix_summary.md     # Ranking-matrix protocol summary
└── requirements.txt        # Exact dependency versions used for the results
```

## Requirements

- Python >= 3.9 (results were produced with Python 3.13.12, CPython, macOS
  arm64 on Apple Silicon hardware)
- numpy 2.4.6, scipy 1.18.0, matplotlib 3.11.1, pillow 12.3.0
  (see `requirements.txt` for the exact reproducibility lock file)
- Optional: pytest >= 8.0 (the test suite also runs standalone without pytest)

## Reproduction

### 1. Install dependencies

```bash
python -m pip install -r requirements.txt
```

### 2. Run the verification suite (optional but recommended)

```bash
python tests/test_verification.py
# or: pytest tests/
```

Every assertion pins one number that appears in the manuscript; a passing run
confirms the environment reproduces the published values.

### 3. Main analysis pipeline

```bash
cd src
python run_all.py
```

This runs the full pipeline with the default seed (20260816) and writes five
figures at 300 dpi, corresponding to the manuscript figures:

| Output file | Manuscript content |
|---|---|
| `Figure1_r0_temperature.png` | R0 as a function of temperature |
| `Figure2_prcc_sensitivity.png` | LHS-PRCC sensitivity bars |
| `Figure3_scenario_peaks.png` | Epidemic peaks across scenarios |
| `Figure4_incidence_timeseries.png` | Incidence time series |
| `Figure5_prcc_designs.png` | PRCC comparison across input-range designs |

### 4. Dynamic-output and ranking-matrix analyses

Run each script from inside `new_analysis/` (the scripts import the model
from `../src`):

```bash
cd new_analysis
python make_static_reference.py    # static R0 PRCC reference
python run_dynamic_prcc.py         # PRCC of dynamic epidemic outputs
python run_scaling.py              # Nv/Nh scaling of the headline peak
python run_clamp.py                # effect of the b(T) clamp
python run_equalwidth.py           # +/-20 % equal-width range comparison
python run_equalwidth50.py         # 50-point equal-width robustness check
python run_static_equalwidth.py    # static R0 under equal-width ranges
python run_constant.py             # constant-parameter PRCC comparison
python run_rank_matrix.py          # ranking-matrix protocol (aggregates the above)
python make_summary.py             # regenerates SUMMARY.md
```

Each script writes its JSON result next to itself; the committed `*.json`
files are the archived outputs used in the manuscript. `SUMMARY.md` and
`rank_matrix_summary.md` give the corresponding headline numbers, rank
comparisons, and Spearman correlations between designs.

### 5. Provenance stamping (optional)

`stamp_src_hashes.py` records the SHA-256 hashes of the source files and a
numeric re-check of the baseline run into the `meta` fields of the result
JSON files. It never modifies any result value.

## Environment

- OS / hardware: macOS (Apple Silicon, arm64)
- Python: 3.13.12 (CPython)
- Locked dependency versions: see `requirements.txt`
- Random seed: 20260816 (fixed; all stochastic steps are seeded)

## Citation

If you use this code, please cite the accompanying manuscript (in
preparation). This repository contains the complete code and archived result
files needed to reproduce every number reported in that study.

## Copyright

Copyright (c) 2026 Tian-Yu Wang. All rights reserved.
