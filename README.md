# Platelet-fibrin transport

Private research distribution of the mechanical model, continuous reflecting
tracer-transport engine, numerical verification, and companion minipaper.
Prepared on 14 October 2026. This is a descriptive research pilot, not a
biologically validated clot predictor.

## Read the compiled documents

- [Minipaper (PDF)](paper/main2vPT.pdf)
- [Numerical supplement (PDF)](paper/validation_supplement.pdf)
- [Five-minute temporal verification guide](docs/temporal_presentation.md)

The PDFs are deliberately versioned alongside their LaTeX sources. Rebuild
them after editing the text or figures; do not treat an old PDF as the current
source. No external Overleaf ZIP is maintained by this repository.

## Two different programs: do not conflate them

1. `src/program19_core/` is the unchanged Program19.3 mechanical code: independent
   two-dimensional layers, saved history, checkpoints, and continuation.
2. `src/continuous_transport/` contains the C++ engine used for the paper's
   **frozen-state, continuous-solid, reflecting** transport and connectivity tools.

The legacy GUI also offers raster-mask, rejected-step and evolving-geometry
transport. Those modes are **not** the current paper's reflecting transport.
See [method and scope](docs/model_and_scope.md).

See the [distribution checks and manuscript path map](docs/repository_checks_20261014.md)
for the checks performed on this snapshot and the locations of its companion data.

## Quick start

Python 3.14.3 with the versions in `requirements-lock.txt` was used for the local
checks. Python 3.11+ may work but is not claimed to give bitwise-identical output.

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip install -e . --no-deps
python scripts/check_repository.py
python -m unittest discover -s tests/mechanics -v
```

Build the current transport engines with a C++17 compiler, then run a tiny
synthetic example (not a new clot simulation):

```sh
python scripts/build_transport.py --compiler g++
python scripts/smoke_transport.py
```

Reproduce the existing cohort's outcome figures from the included per-particle
event records, without rerunning mechanics or transport:

```sh
python scripts/reproduce_paper_figures.py
```

Generated figures go into `build/reproduced/`, not over the published figures.
The historic plotting routine is reused unchanged. Its generated table predates
later editorial row separators; the canonical typeset table is under `paper/`.

Compile both documents (requires a LaTeX installation with `latexmk`, `pdflatex`
and `bibtex`):

```sh
python scripts/build_paper.py
```

This refreshes the two tracked PDFs in `paper/`; temporary TeX output stays in
`paper/build/`. Inspect them before committing.

## Mechanics

```sh
python src/program19_3.py --help
python src/program19_3.py
```

The second command opens the GUI. Use **8 layer processes** as the working
default on the originating PC. Hardware and memory requirements vary. Running
the full mechanical cohort is an expensive, separate operation, not a prerequisite
for viewing the paper or reproducing the included outcome plots.

## What is included

- Program19.3 source and its existing unit tests.
- Exact archived continuous-transport C++ source, including the separate finer
  Brownian-increment version used in the temporal study.
- 29,952 production event histories over 117 unique geometries, with paired
  contrasts and geometry summaries (nine mechanical seed groups).
- The 73,728 temporal-verification event histories and original/revised criteria.
- LaTeX, a single shared bibliography, required figures, and compiled PDFs.
- Import hashes and source provenance in `data/file_manifest.json`.

## What is not included

The large mechanical trajectories, checkpoints, and complete primitive scenes
have **not yet been deposited**. There is no invented download URL or DOI.
The full mechanical-to-transport workflow is therefore not yet a one-command
reproduction from this repository alone. See [reproducibility](docs/reproducibility.md)
and [data inventory](data/README.md).

`archive/research_drivers/` preserves original research scripts, including their
historical paths. It is evidence of the workflow, not the portable entry point.
Do not run these scripts unmodified against someone else's filesystem.

## Publication and licensing status

Keep this repository private until author approval, third-party licensing, data
archiving, and metadata have been reviewed. No blanket open-source license has
been invented. See [license status](LICENSE_STATUS.md) and
[third-party attribution](THIRD_PARTY_NOTICES.md).

Repository contents and PDFs do not establish complete temporal/spatial or
biological validation. The temporal acceptance uses a **post-analysis revision**
from ±2 to ±2.5 percentage points, retaining ±3% for relative RMST. The original
decision is preserved, not relabeled.
