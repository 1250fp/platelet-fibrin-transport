# Distribution checks - 14 October 2026

This record concerns packaging and regression checks, not a new scientific
validation or a new mechanical/transport cohort.

| Check | Result |
|---|---|
| Imported file integrity | 382 SHA-256 hashes matched the assembly manifest |
| Production event records | 117 unique geometries, 29,952 histories; S(100) and RMST(100) recomputed from event times |
| Temporal verification records | Nine CSVs, 8,192 histories each; 73,728 histories kept separate from production |
| Existing Python tests | 37 passed; mechanics and legacy raster-transport tests |
| C++ builds | Four current transport/connectivity/clearance executables compiled |
| Synthetic reflecting-transport smoke test | Passed deterministic solid reflection, free motion, axial exit and 16 short Brownian histories |
| Outcome figure reconstruction | Completed from included event records, without rerunning the cohort |
| LaTeX builds | Main paper: 14 pages; supplement: 13 pages; no unresolved references or overfull-box warnings in final logs |
| PDF visual check | All pages reviewed in contact sheets; parameter table and revised temporal criteria inspected at higher resolution |
| Secret-pattern check | No GitHub-token, OpenAI-key or private-key markers found in the distribution |

The secret-pattern check is a basic precaution, not a comprehensive security
audit. Original local paths remain in archived scripts and provenance records.
Compiled executables, build logs, Python caches and the large raw mechanical
archive are excluded from Git.

## Environment

- Windows; Python 3.14.3 for scientific tests and plotting.
- NumPy 2.4.3, SciPy 1.17.1, Matplotlib 3.10.8, Pillow 12.1.1.
- GCC via Strawberry Perl's C++ toolchain; MiKTeX/latexmk for both documents.
- A separate bundled Python was used for standard-library integrity checks.

No numerical parameters or simulation source algorithms were changed for
packaging. The paper copy uses one shared bibliography. Markdown links use
repository-relative paths. Original source directories and the Overleaf ZIP
were not modified.

## Paths printed in the historical manuscripts

The manuscripts retain their original companion-data labels. In this repository,
the corresponding files have been organized as follows (paths start at the
repository root):

| Printed/original companion path | Repository location |
|---|---|
| `data/homogeneous_pilot/` | `data/processed/homogeneous_pilot/` |
| `data/validation/` | `validation/historical_data/` |
| `data/temporal_confirmation_20261009/` | `validation/temporal_20261009/` and `validation/temporal_20261012_revision/` |
| `figures/validation/` | `paper/figures/validation/` |

The small `.tex` table inputs remain under `paper/data/` so both papers compile
without an external data directory. They are not the complete scientific data
archive. See [reproducibility](reproducibility.md) for the remaining external
inputs needed to rerun the full workflow.
