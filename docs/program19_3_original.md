# Program19.3

Program19.3 extends the audited Program19.2 mechanics with three non-destructive
workflows:

1. continue a completed mechanical stack for additional physical time;
2. add independent mechanical layers with new, non-colliding seeds; and
3. release Brownian tracers at the stack midpoint and measure first passage to
   either absorbing z face while the saved clot geometry evolves.

It does **not** generate autoregressive (AR) layers. Every layer is an
independent two-dimensional mechanical realization and the ordered collection
is interpreted as a spatial z stack.

## Parallel layer execution

Mechanical layers can run concurrently because they have different seeds and
do not exchange forces. This applies to new stacks, added layers, and temporal
continuations. In the GUI, **Parallel layer processes** selects the maximum
number of simultaneous layers; continuation and extension dialogs also ask for
the desired value. A value of `1` preserves serial execution. Values greater
than one use separate Python processes rather than threads, allowing CPU-bound
mechanics to use multiple CPU cores despite Python's global interpreter lock.
The effective count is capped by the layer count and the logical CPUs visible
to Python.

Each process limits numerical-library worker threads to one, preventing nested
BLAS/OpenMP oversubscription. Start with 2--4 processes; increasing the value
also increases peak RAM use and may stop improving speed once memory bandwidth
or CPU cooling becomes limiting. Parallel scheduling does not change the seed
rule, equations, saved layer order, or the absence of interlayer forces.

## Starting the GUI

Double-click `run_program19_3.bat`, or run:

```powershell
py -3.14 program19_3.py
```

The **Previous runs and evolving transport** tab scans these libraries by
default:

- `program19_2/results`
- `program19_3/results`

Additional libraries can be appended with **Browse...**. Select a row and use:

- **Inspect parameters**: view the principal mechanical and transport settings.
- **Continue mechanical evolution...**: evolve every saved layer from its final
  state for an additional time without regenerating its initial geometry.
- **Extend stack...**: reuse the existing layers and generate only the missing
  independent layers with new seeds.
- **Evolving transport...**: release tracers at the z midpoint and compare the
  evolving geometry with the same geometry frozen at release time.
- **Frozen transport...**: select any saved mechanical time and reuse that fixed
  stack for a post-hoc transport horizon independent of the mechanical duration.

Every operation creates a timestamped child directory. Existing simulations
are never overwritten.

## Mechanical continuation and reproducibility

Continuation restores node and platelet positions, topology, WLC rest data,
cohesive links and their birth/reference lengths, filopodia, occupied-node
state, simulation time, and step count.

New Program19.3 runs save `checkpoint.npz`, including the NumPy random-number
generator state. Continuations from these checkpoints therefore preserve the
stochastic sequence exactly. Program19.2 did not save that RNG state;
continuing a legacy run consequently starts a documented deterministic legacy
branch from the saved physical state. It is reproducible, but it is not the
unrecoverable random sequence that the old process would have drawn.

The child run contains:

- `continuation_manifest.json`, with parent, time interval, and RNG policy;
- one folder per continued layer;
- `trajectory.npz` and `intermediate_masks.npz` when those outputs are enabled;
- `mechanical_spatial_stack.npz`, pore/accessibility summaries, figures, and
  status files at the new final time.

Temporal continuation currently requires a completed stack whose layers share
the same final time. A cancelled partial continuation remains saved, but is not
automatically resumed by this version.

## Extending the spatial stack

Stack extension does not duplicate the original large layer files. The new
result stores a `stack_manifest.json` that references reused layers and newly
generated layers, plus a combined final `mechanical_spatial_stack.npz`.
Consequently, keep the referenced parent and `_components` directories with the
extended result. Seeds are checked against the configured simulation
libraries before generation.

To combine both operations, first extend the original stack and then continue
the enlarged stack in time. Extending a temporal-continuation child directly is
intentionally rejected because its parent layers and continued layers represent
different time intervals.

## Evolving-clot transport

The transport experiment has one deliberately narrow interpretation:

- all tracers start in accessible voxels near the stack midpoint;
- `z=0` and `z=L` are both absorbing escape faces;
- lateral boundaries are reflecting;
- attempted moves into the solid phase are rejected;
- tracers do not interact with one another and do not alter the clot;
- there is no inlet-to-outlet injection experiment in this analysis.

For the evolving mode, the solid and analyzed-region masks are piecewise
constant between saved mechanical frames (normally separated by the configured
output interval). If a changing boundary covers a tracer center, the center is
projected to the nearest accessible voxel and the event is counted as a
`moving_boundary_relocation`; it is not interpreted as biological capture or
aggregation.

An evolving experiment requires mechanical history through

```text
release time + observation time.
```

Continue the mechanical run first if the saved history is too short. The
optional frozen comparison uses the same release geometry, particle initial
conditions, seed, and Brownian increments, but does not update the mask. This
isolates the effect of geometry evolution from Monte Carlo noise.

Transport output contains:

- `analysis_manifest.json`;
- `evolving_summary.json` and `evolving_first_passage.npz`;
- corresponding `frozen_*` files when requested;
- `comparison_summary.json` and `survival_comparison.png`.

Reported quantities include escape fractions through each z face, censored
fraction, restricted mean first-passage time, survival curve, rejected solid
and lateral moves, and moving-boundary relocations. The saved trajectory subset
is for visualization; summary statistics use every simulated tracer.

## Command-line examples

List reusable runs:

```powershell
py -3.14 program19_3.py --catalog-root "E:\path\to\results"
```

Continue mechanics by two seconds:

```powershell
py -3.14 program19_3.py --continue-run "E:\path\to\run" `
  --additional-time 2 --workers 4 `
  --output-root "E:\path\to\program19_3_results"
```

Extend a five-layer stack to ten layers:

```powershell
py -3.14 program19_3.py --extend-stack "E:\path\to\run" `
  --target-layers 10 --first-new-seed 12035 `
  --workers 4 `
  --output-root "E:\path\to\program19_3_results"
```

Run central-release transport from mechanical time 2 s for a 1 s observation:

```powershell
py -3.14 program19_3.py --evolving-transport-run "E:\path\to\run" `
  --release-time 2 --observation-time 1 --particle-radius 0.05 `
  --diffusion 4.0 --particles 500 --transport-seed 2035 `
  --output-root "E:\path\to\transport_analyses"
```

Use `--no-frozen-comparison` when only the evolving geometry is needed.

Run frozen central-release transport from the saved mechanical frame at 2 s and
observe it for 20 s, without continuing or regenerating the mechanical stack:

```powershell
py -3.14 program19_3.py --frozen-transport-run "E:\path\to\run" `
  --release-time 2 --observation-time 20 --particle-radius 0.05 `
  --diffusion 4.0 --particles 500 --transport-seed 2035 `
  --output-root "E:\path\to\transport_analyses"
```

Generate ten independent layers with four parallel processes:

```powershell
py -3.14 program19_3.py --headless --layers 10 --workers 4
```

## Verification

Run the test suite with the Python installation that contains NumPy, SciPy,
and Matplotlib:

```powershell
py -3.14 -m unittest discover -s tests -v
```

The suite checks the audited WLC and bending behavior, reciprocal platelet
interactions, cohesion rules, corrected diffusion, profile isolation,
checkpoint/RNG restoration, completed-run continuation, spatial extension, and
the time-dependent first-passage implementation.
