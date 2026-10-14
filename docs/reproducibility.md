# Reproducibility levels

## Available now

1. Compile the paper and numerical supplement from this checkout alone, with TeX.
2. Recompute S(100), RMST(100), exit counts and cohort outcome figures from the
   included per-particle CSVs. Preserve 117 unique geometries / 29,952 histories.
3. Inspect the numerical-study outcomes, paired intervals, original protocol and
   explicitly revised tolerance. Original JSONs are retained without editing.
4. Build the archived C++ sources and run the small synthetic smoke example.
5. Run the existing Program19.3 unit tests, including mechanics and legacy transport.

The provided smoke example is an installation/regression check, not a substitute
for the full numerical-verification suite. The checked-in plots and original
validation reports document the completed studies.

## Requires the external mechanical archive

Recreating all continuous scenes, continuing saved mechanics, re-certifying
connectivity or repeating the complete cohort requires the original trajectories,
checkpoint states, support rasters and referenced parent stack components. These
large objects are intentionally not included in this first distribution.

Historical runners under `archive/research_drivers/` preserve source logic and
original private paths. They require path remapping and external inputs. Only
the scripts listed in the root README are portable entry points in this snapshot.
Independent figures for the mechanical evolution / 50-layer solid rendering are
included as finished assets; regenerating those requires the external histories.

## PDF policy

Both `paper/main2vPT.pdf` and `paper/validation_supplement.pdf` are tracked.
Run `python scripts/build_paper.py` after source edits. Commit the changed sources,
required figures/tables and rebuilt PDF together. Build products other than the
two final PDFs are ignored. PDF bytes can differ due to TeX version, metadata and
build date; numerical tables must not change merely because of recompilation.

## Version policy

This is the initial private preparation, not a published release. Assign a release
tag and persistent data identifier only after checking author metadata, licensing,
and the deposit. Keep future stereology/density revisions distinct from this cohort.
