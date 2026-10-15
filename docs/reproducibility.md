# Reproducibility and verification record

## Available from this checkout

1. Compile the minipaper and numerical supplement with a TeX installation.
2. Recompute S(100), RMST(100), exit counts and outcome figures from the included
   117 unique production geometries / 29,952 per-particle histories.
3. Inspect numerical-study outcomes, paired intervals, the original protocol and
   the explicitly revised temporal tolerance.
4. Build the continuous C++ tools and run the small synthetic smoke example.
5. Run Program19.3 unit tests for mechanics and legacy raster transport.

The smoke example is an installation/regression check, not a substitute for the
numerical-verification studies. Original scientific values are retained; personal
paths in protocols and provenance were normalized for distribution.

## Requires external inputs

Recreating continuous scenes, continuing saved mechanics, re-certifying
connectivity or rerunning the cohort requires original trajectories, checkpoints,
support rasters and referenced parent stack components. These are not deposited
in this repository. The finished mechanical-evolution and 50-layer solid figures
also require those histories to regenerate.

The historical runners in `archive/research_drivers/` require path remapping and
external inputs. Only commands documented in the root README are portable entry
points. [External-data identifiers and hash semantics](../data/README.md) are
documented separately. The full cohort is not yet a one-command reproduction.

## Distribution checks: 14 October 2026

These are the recorded checks on the initial snapshot, not new scientific
validation or a new mechanical/transport cohort.

| Check | Recorded result |
|---|---|
| Imported integrity | 382 SHA-256 hashes matched the initial assembly manifest |
| Production events | 117 unique geometries, 29,952 histories; S(100) and RMST(100) recomputed from event times |
| Temporal events | Nine CSVs, 8,192 histories each; 73,728 histories kept separate from production |
| Python tests | 37 passed; mechanics and legacy raster-transport tests |
| C++ builds | Four transport/connectivity/clearance executables compiled |
| Synthetic smoke test | Deterministic reflection, free motion, axial exit and 16 short Brownian histories passed |
| Outcome figures | Reconstructed from included events without rerunning the cohort |
| LaTeX | Main paper 14 pages, supplement 13 pages; no unresolved references or overfull boxes in final logs |
| PDF review | All pages reviewed; parameter table and temporal criteria inspected at higher resolution |
| Basic secret-pattern scan | No GitHub-token, OpenAI-key or private-key markers detected |

Environment: Windows, Python 3.14.3, NumPy 2.4.3, SciPy 1.17.1, Matplotlib 3.10.8,
Pillow 12.1.1; GCC from Strawberry Perl and MiKTeX/latexmk. The requirements file
pins direct scientific dependencies, not a complete cross-platform environment.
A separate bundled Python performed standard-library integrity checks. These
results do not assert bitwise reproducibility with other libraries or compilers.

## Distribution cleanup: 15 October 2026

Personal presentation notes were retained outside the repository. Duplicate
historical protocols were consolidated, documentation merged, local source roots
replaced by logical labels, and missing original temporal-study figures restored.
Original and revised temporal criteria remain distinct. No scientific simulation
from the production or validation cohorts was rerun; no numerical data,
mechanics/transport algorithms, paper sources or
compiled manuscript PDFs were revised by this cleanup.

The manifest now records original and distributed hashes separately. Its retired
records and aliases document removals and deduplication. `scripts/check_repository.py`
checks distributed hashes, cohort outcomes, relative Markdown links and accidental
personal-path leakage. This is a packaging safeguard, not a comprehensive security
audit. Earlier Git commits retain their original content; history was not rewritten.

Checks after cleanup used an isolated Python 3.12.14 environment with the pinned
direct scientific dependencies. All 37 unit tests passed, the 16-history synthetic
C++ smoke example passed, and outcome figures were regenerated under ignored
`build/` without replacing published assets. All 374 distributed manifest hashes
and relative Markdown links passed. A comparison against the pre-cleanup snapshot
confirmed 269 byte-identical CSVs, one inventory CSV with path-only changes, 35
JSON files with unchanged content apart from path labels (including deduplicated
copies), 38 unchanged source/test files, and 19 unchanged manuscript assets.
The 13 historical Python drivers retain the same syntax trees apart from source
path strings. Tracked content decreased from 404 files / 48.95 MiB to 394 files /
44.46 MiB before Git compression. These checks do not broaden scientific validation.

## Companion paths printed in the manuscripts

The manuscripts retain their original companion-data labels. Locations below
start at the repository root:

| Printed/original label | Repository location |
|---|---|
| `data/homogeneous_pilot/` | `data/processed/homogeneous_pilot/` |
| `data/validation/` | `validation/historical_data/` |
| `data/temporal_confirmation_20261009/` | `validation/temporal_20261009/` and `validation/temporal_20261012_revision/` |
| `figures/validation/` | `paper/figures/validation/` |

Small `.tex` inputs remain under `paper/data/` for compilation; they are not the
complete scientific archive.

## PDF and version policy

Both `paper/main2vPT.pdf` and `paper/validation_supplement.pdf` are tracked.
After changing manuscript sources or figures, run `python scripts/build_paper.py`,
inspect the output, and commit sources and PDFs together. Temporary TeX files
are ignored. Different TeX metadata can change PDF bytes without changing results.
The repository does not update an external Overleaf ZIP.

This is a public research snapshot, not a fully archived journal/data release.
Licensing and author metadata remain under review; no license or DOI is inferred
from public visibility. Keep future stereology/density revisions distinct from
the existing cohort, and retain the original +/-2 pp temporal decision alongside
the post-analysis +/-2.5 pp revision.
