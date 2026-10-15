# Included data and provenance

`processed/homogeneous_pilot/` contains production event records, per-state
connectivity/escape classifications, geometry summaries, paired contrasts,
survival curves and source protocols. Its 117 canonical state folders contain
29,952 histories. The 144 nominal force/time labels include shared initial
geometries; do not count those histories repeatedly.

Temporal verification is separate under `validation/temporal_20261009/`, with
nine event CSVs and 73,728 histories. The original +/-2 pp criterion remains
separate from its post-analysis revision under `validation/temporal_20261012_revision/`.
Neither these histories nor their conditional intervals are extra mechanical
replicates in the production cohort.

## External-data labels

Paths beginning with `external_data/` are logical source identifiers, not download
links or files supplied by this checkout. They replace personal filesystem roots:

| Logical root | Meaning |
|---|---|
| `external_data/research_workspace/` | Historical preparation scripts, outputs and analyses |
| `external_data/program19_3/` | Original mechanics source, histories and checkpoints |
| `external_data/manuscript/` | Original manuscript work area and prior snapshots |
| `external_data/fibrin_project/` | Other historical project inputs |

For included files, use repository-relative paths. To rerun a historical batch,
obtain the external inputs and explicitly remap its labels; the prefix does not
make the archive available. Both `external_data/` and `data/raw/` are ignored by
Git to avoid accidental uploads. No raw-data deposit or DOI is claimed.

## Manifest and deduplication

`file_manifest.json` distinguishes original imported bytes (`source_sha256` and
`source_bytes`) from current distributed bytes (`sha256` and `bytes`).
`assembly_sha256` and `assembly_bytes` preserve the initial repository snapshot
where available. Original source sizes are null when they were not recorded or
could not be reverified. Hashes inside
historical protocols still refer to original inputs, not automatically to their
path-normalized copies. Path labels change; scientific values and original hashes
do not. `retired_records` preserves the identity of removed personal notes, while
`aliases` records former duplicate locations and renamed documentation.

The canonical copies of replica 5, replica 6, replicas 1-2 and replicas 3-4 source
protocols and execution summaries are under
`processed/homogeneous_pilot/source_protocols/historical_batches/`. The duplicate
top-level copies were removed after byte-for-byte verification. New-replica/15 s
and previous-six-replica consolidation protocols remain separate.

Event CSVs, numerical results, original/revised criteria and scientific figures
are retained. The small `.tex` table inputs under `paper/data/` are intentionally
available to compile the documents independently of analysis scripts.
