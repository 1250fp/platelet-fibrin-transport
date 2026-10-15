# Platelet-fibrin transport

Research code, processed simulation data, numerical verification, and companion
minipaper on tracer retention in mechanically reorganized platelet-fibrin stacks.
This is a descriptive pilot, not a biologically validated clot predictor.

- [Minipaper (PDF)](paper/main2vPT.pdf)
- [Numerical supplement (PDF)](paper/validation_supplement.pdf)
- [Model and scientific scope](docs/model_and_scope.md)
- [Mechanical simulation guide](docs/mechanics.md)
- [Reproducibility and verification record](docs/reproducibility.md)
- [Data inventory and provenance](data/README.md)

## Two distinct workflows

`src/program19_core/` generates independent 2D mechanical layers, histories and
checkpoints. `src/continuous_transport/` contains the C++ engines used for the
paper's **frozen-state, continuous-solid, reflecting** transport and connectivity.
The GUI's legacy raster-mask and evolving-geometry transport is **not** the
paper's reflecting transport. See [model and scope](docs/model_and_scope.md).

## Install and check

Run from the repository root. Python 3.14.3 was used for the recorded scientific
checks. Python 3.11+ is supported by the package metadata but bitwise-identical
output across environments is not claimed.

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip install -e . --no-deps
python scripts/check_repository.py
python -m unittest discover -s tests/mechanics -v
```

## Use

```sh
# Inspect options or open the mechanical GUI
python src/program19_3.py --help
python src/program19_3.py

# Build the transport tools (requires a C++17 compiler) and check a tiny example
python scripts/build_transport.py --compiler g++
python scripts/smoke_transport.py

# Recreate outcome plots from included events; no simulation is rerun
python scripts/reproduce_paper_figures.py

# Compile both PDFs (requires latexmk, pdflatex and bibtex)
python scripts/build_paper.py
```

Generated outcome figures go to `build/reproduced/`. Paper compilation updates
the two tracked PDFs in `paper/`; inspect them before committing. Eight layer
processes is the working choice on the originating PC, not a universal hardware
requirement or the code's default. See the [mechanics guide](docs/mechanics.md).

## Included data and limits

The distribution includes 29,952 production histories over 117 unique geometries
(nine mechanical seed groups), plus 73,728 separate temporal-verification histories.
Full mechanical trajectories, checkpoints and continuous primitive scenes have
**not yet been deposited**. This is not a one-command reproduction of the full
mechanical-to-transport cohort. [Available and missing inputs](docs/reproducibility.md)
are distinguished explicitly; no unavailable download or DOI is implied.

`archive/research_drivers/` retains historical workflow logic, with personal roots
replaced by external-data labels. These are not portable batch entry points;
the documented scripts above are. No external Overleaf ZIP is maintained here.

## Scientific and licensing status

Temporal acceptance uses a **post-analysis revision** from +/-2 to +/-2.5
percentage points for survival, retaining +/-3% for relative RMST. Original
results and the original inconclusive decision are preserved separately from
the revised criterion. Neither this test nor the repository establishes complete
spatial, all-state temporal, or biological validation.

The repository is publicly visible. Ownership, collaborator approval and
third-party licensing remain under review; no blanket open-source license is
assigned. See [license status](LICENSE_STATUS.md) and
[third-party attribution](THIRD_PARTY_NOTICES.md).
