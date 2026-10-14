# Included data

`processed/homogeneous_pilot/` contains production event records, per-state
connectivity/escape classifications, geometry summaries, paired contrasts,
survival curves and source protocols. They are retained with their original
numerical content. The 117 canonical state folders contain 29,952 histories.

`file_manifest.json` records source paths and hashes at distribution assembly.
Historical absolute paths in provenance describe where data came from; portable
scripts resolve included data relative to this repository.

The temporal verification is separate under `validation/temporal_20261009/`.
Its nine merged CSVs contain 73,728 histories. These must not be added to the
production cohort. The post-analysis margin revision has its own directory.

`data/raw/` is ignored and reserved for a future external data archive. No archive
has yet been uploaded and no DOI or download URL is claimed. Original full
mechanical trajectories, checkpoints and primitive scenes remain in the author's
local work area. Original source data are not deleted by packaging.
