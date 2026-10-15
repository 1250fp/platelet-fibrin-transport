# Mechanical simulation guide

Program19.3 constructs independent two-dimensional platelet-fibrin layers and
organizes them as a spatial stack. Layers have separate seeds and do not exchange
mechanical forces. They are not autoregressive copies or a coupled 3D mechanics
simulation. Run the commands below from the repository root after installation.

## Start and configure

```sh
python src/program19_3.py --help
python src/program19_3.py
```

The second command opens the GUI. Configuration distinguishes equation,
parameter, numerical and geometry profiles. A profile name alone is not a full
record of a run: preserve its `config.json` and `manifest.json`.
`configs/legacy_program19_example.json` is a legacy example, not the current
paper's transport protocol.

For a new mechanical run, select duration, output interval, layer count, seeds,
platelet count and physical parameters. This illustrative command generates
eight independent layers for one second; it does not reproduce the paper cohort:

```sh
python src/program19_3.py --headless --layers 8 --workers 8 --duration 1 --seed 4000 --output-root results/mechanics
```

## Parallelism and paired runs

`--workers` sets the number of independent layer processes, not threads working
on one layer. Eight processes is the chosen operating setting on the originating
PC; the code default remains one. The effective count is capped by layer count
and available logical CPUs. More processes increase memory use. When available,
`threadpoolctl` limits numerical-library threads inside each process.

A layer uses the base seed plus its index times the configured seed stride.
The sweep option `--f0-sweep` shares the same seed block across force conditions
within a replicate; `--replicates` assigns distinct layer-seed blocks to subsequent
replicates. This pairing is distinct from adding more layers to one stack.

## State and output

Each layer directory records its configuration, initial and final states,
diagnostics and checkpoint. With trajectory saving enabled, `trajectory.npz`
contains times and positions at the output interval; `trajectory_detail=topology`
also retains time-dependent cohesive links and filopodia. The checkpoint stores
the mechanical state and random-number generator state for continuation.

The integrator may use adaptive internal substeps. The output interval is a
storage interval, not the Brownian or force-integration step. Periodic checkpoints
replace the latest checkpoint; they are not a separate restart file for every
saved trajectory frame.

At stack level, the run records an index of layers, a manifest, configuration,
assembled masks and summaries. Keep the manifest and parent references together
with the layer files.

## Find and reuse existing runs

The GUI's previous-runs tab can scan additional result folders. From the CLI:

```sh
python src/program19_3.py --catalog-root results
```

Replace `results/example_run` below with an actual completed run. Continue its
existing layers for two additional seconds:

```sh
python src/program19_3.py --continue-run results/example_run --additional-time 2 --workers 8 --output-root results/continued
```

Or extend a smaller stack to 50 layers with an unused seed block:

```sh
python src/program19_3.py --extend-stack results/example_run --target-layers 50 --first-new-seed 12000 --workers 8 --output-root results/extended
```

These operations create child directories rather than overwriting the source.
Extension references reused layers instead of duplicating their large files:
retain the parent and component directories. Extend the original stack before
continuing it in time; direct extension of a temporal-continuation child is
rejected. A cancelled partial continuation is not automatically resumable by
this version's completed-stack continuation command.

New checkpoints restore the saved random-generator state. Legacy inputs without
it use a documented deterministic fallback seed; that is reproducible, but it
does not recover the interrupted run's missing random sequence.

## Module map

| Responsibility | Modules in `src/program19_core/` |
|---|---|
| Interface and settings | `cli`, `gui`, `config`, `profiles` |
| Initial geometry and state | `geometry`, `codeocean_geometry`, `state` |
| Forces and integration | `wlc`, `bending`, `interactions`, `cohesion`, `integrator`, `spatial` |
| Execution and storage | `mechanical_runner`, `parallel`, `persistence`, `reporting` |
| Reuse and experiments | `catalog`, `continuation`, `stack_extension`, `experiments` |
| Legacy raster analysis | `layers`, `transport`, `transport_steps`, `evolving_transport`, `through_clot`, `pilot` |
| Package/compatibility entry points | `__init__`, `runner` |

## Transport is a separate methodological choice

The GUI and CLI expose legacy raster-based first-passage and evolving-geometry
analysis. Their masks, rejected moves and moving-boundary handling must not be
conflated with the paper's continuous reflecting C++ engine. The basic legacy
rasterizer uses backbone pairs; the paper's preparation separately includes
cohesive segments as obstacles. Saved topology is therefore important when
reconstructing the paper geometry.

The paper evaluates frozen mechanical states. The full preparation requires the
external mechanical archive; see [reproducibility](reproducibility.md). No new
simulation is required to reproduce its included outcome plots.

## Tests

```sh
python -m unittest discover -s tests/mechanics -v
```

These tests cover mechanics, profile separation, saved state, continuation,
stack extension and legacy raster transport. They do not by themselves certify
the separate continuous-transport numerical studies or biological parameters.
