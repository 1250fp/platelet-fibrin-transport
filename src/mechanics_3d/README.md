# fibrin-clot-sim-cuda

3D CUDA simulator of platelet-driven contraction of fibrin clots: a fibrin network (worm-like-chain fibres with torsion, dynamic cross-linking) plus platelets whose filopodia pull on fibres and other platelets.

This is a **modified version** of the model by Christian Michael, Francesco Pancaldi and co-authors (code copyright (c) 2022 Zhiliang Xu, MIT License, see `LICENSE`). Original code: Code Ocean capsule 2690247. If you use it, cite the original paper:

> C. Michael, F. Pancaldi, S. Britton, O. V. Kim, A. D. Peshkova, K. Vo, Z. Xu, R. I. Litvinov, J. W. Weisel, M. Alber. *Combined computational modeling and experimental study of the biomechanical mechanisms of platelet-driven contraction of fibrin clots.* Communications Biology 6, 869 (2023). https://doi.org/10.1038/s42003-023-05240-z

Modifications by Paolo Teppa.

## Changes from the original

Only 5 source files differ from Michael's code. Formulas, parameters, integrator and random-number code are unchanged.

1. **Force-bound fix** (`src/functor_misc.h`, lines 34 and 71): `maxNodeCount(maxNodeCount)` changed to `maxNodeCount(_maxNodeCount)`. The original initialised the bound from itself, an uninitialized value, in the functors that add torsion, platelet arm-node, platelet field-node and volume-exclusion forces to the nodes. **This can change results** relative to the original.
2. **Params export fix** (`src/src_fibrin/functor_calc_strain_params.h`): edge strain `(L-L0)/L0` and alignment `|dz/L|` are looked up in the correct neighbour slot. Invalid entries are written as NaN. In the original these output columns were mostly garbage. Output only, no effect on the dynamics.
3. **Optional stop time** (`src/System.cu`): the environment variable `MICHAEL_FINAL_TIME_SECONDS` sets the simulated end time. Without it, runs stop at the original 2700 s (45 min).
4. **Faster output** (`src/Storage.cu`): one bulk device-to-host copy per snapshot, and `'\n'` instead of `std::endl`. Written values and formatting are unchanged.
5. **Thrust compatibility** (`src/SystemStructures.h`): `thrust::get<0>(x)` and an added include, needed by CUDA 12's Thrust.
6. **Build**: a single `CMakeLists.txt` that builds Release (`-O3`), C++14, CUDA 12, without fast-math. The original was a CUDA 10 debug build (`-g -G`). Expect statistical agreement with the original, not identical trajectories.

Removed from the original tree: per-folder Makefiles and CMake files, editor swap files, sample XMLs, unused glm headers, and pugixml docs, scripts and prebuilt library. The original reset or resume mode (`bend-model reset ...`) is still in `main.cu`, but it does not reliably restore every state field, so start runs from an XML.

**Known bug (not fixed):** the Params file header writes `minZ`/`maxZ` from `minX`/`maxX` (`src/Storage.cu` lines 278–279). The State-file header is correct.

## Build (Linux or WSL2, NVIDIA GPU)

You need the CUDA toolkit 12.x (12.8 or newer for RTX 50-series), CMake ≥ 3.24 and g++.

```bash
cmake -S . -B build            # default GPU arch 120 (RTX 50-series)
cmake --build build -j4         # each nvcc job needs ~1-2 GB RAM; plain -j can run out of memory
# other GPUs: cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=89   (RTX 40; 86 = RTX 30)
```

This produces `build/bend-model`.

## Run

```bash
./run.sh examples/cylinder_R24_H32_seed1.xml          # 1800 simulated s, dt = 0.001
./run.sh my_clot.xml 600 /path/to/output              # 600 s, custom folder
```

`run.sh` creates the output folder (default `runs/<input>_<T>s/`, which must be empty), copies the XML into it as `input.xml`, then runs the solver there and writes `execution.log` and `exit_code.txt`. To call the solver directly:

```bash
MICHAEL_FINAL_TIME_SECONDS=1800 ./build/bend-model -dt=0.001 input.xml
```

The XML path must be the last argument. `-dt` defaults to 0.001.

## Input

An XML file containing `<settings>` (model parameters), `<nodes>` (fibrin node `x y z`), `<plts>` (platelet centres), `<links>` (node index pairs) and `<fixed>`. The clot shape comes only from the coordinates. The three files in `examples/` are cylinder clots (R = 24 µm, H = 32 µm, seeds 1–3) made with the companion generator `cylinder-xml-generator`, using Michael's "Standard" parameter set.

## Output

The solver writes three folders into the working directory, named after the run settings (`<kind>__domain_..._plt_count_..._filo_count_...`):
- `State__*/State_<t>.sta`: full-state snapshots.
- `Params__*/Param_<t>.sta`: per-snapshot text summaries (time, domain, node, platelet and edge data).
- `Animation__*/`: `.vtk` files for ParaView (`FibrinNetwork_NNNNN.vtk`, `PlateletNNNNN.vtk`, `PlateletConnNNNNN.vtk`).

`tools/make_vtk_series.py` (Python 3, standard library only) writes `.vtk.series` indexes so ParaView shows the real simulated times. Run it after a successful run:

```bash
python3 tools/make_vtk_series.py runs/<run_folder> --expected-final-time 1800
```

## Third-party code

- `src/src_ext_lib/pugixml/`: pugixml 1.9, MIT (c) Arseny Kapoulkine. See its `readme.txt`.
- `src/src_ext_lib/glm/`: GLM 0.9.9.5, MIT (c) G-Truc Creation. See `copying.txt`.
