"""Discovery and concise metadata for Program19.2/19.3 mechanical runs."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class SimulationRecord:
    path: str
    program: str
    created_local: str
    status: str
    completed_layers: int
    requested_layers: int
    mechanical_time_s: float
    f0_nN: float
    base_seed: int
    layer_spacing_um: float
    equation_profile: str
    parameter_profile: str
    numerical_profile: str
    has_time_history: bool
    exact_rng_checkpoint: bool
    parent_path: str | None

    def to_dict(self) -> dict:
        return asdict(self)


def _read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return {}


def _layer_directories(run_directory: Path) -> list[Path]:
    stack_manifest = run_directory / "stack_manifest.json"
    if stack_manifest.is_file():
        data = _read_json(stack_manifest)
        return [Path(value).resolve() for value in data.get("layer_directories", [])]
    directory = run_directory / "mechanical_layers"
    if not directory.is_dir():
        return []
    return sorted(
        path for path in directory.iterdir()
        if path.is_dir() and path.name.startswith("layer_")
    )


def inspect_run(run_directory: str | Path) -> SimulationRecord | None:
    run = Path(run_directory).resolve()
    if "_components" in run.parts:
        return None
    if not (run / "config.json").is_file():
        return None
    if not ((run / "mechanical_spatial_stack.npz").is_file() or _layer_directories(run)):
        return None

    config = _read_json(run / "config.json")
    manifest = _read_json(run / "manifest.json")
    status_data = _read_json(run / "status.json")
    continuation = _read_json(run / "continuation_manifest.json")
    stack_manifest = _read_json(run / "stack_manifest.json")
    layers = _layer_directories(run)
    summaries = [_read_json(path / "layer_summary.json") for path in layers]
    times = [
        float(item.get("simulated_time_s", 0.0))
        for item in summaries if item
    ]
    histories = [path / "intermediate_masks.npz" for path in layers]
    checkpoints = [path / "checkpoint.npz" for path in layers]
    created = str(manifest.get("created_local") or continuation.get("created_local") or "")
    if not created:
        created = datetime.fromtimestamp(run.stat().st_mtime).isoformat(timespec="seconds")

    requested = int(
        status_data.get("requested_mechanical_layers", config.get("mechanical_layers", len(layers)))
    )
    completed = int(
        status_data.get("completed_or_partial_layers_saved", len(layers))
    )
    return SimulationRecord(
        path=str(run),
        program=str(manifest.get("program", "Program19.2 legacy")),
        created_local=created,
        status=str(status_data.get("status", "unknown")),
        completed_layers=completed,
        requested_layers=requested,
        mechanical_time_s=min(times) if times else float(config.get("duration_s", 0.0)),
        f0_nN=float(config.get("filopodia_f0_nN", 0.0)),
        base_seed=int(config.get("seed", 0)),
        layer_spacing_um=float(config.get("layer_spacing_um", 0.0)),
        equation_profile=str(config.get("equation_profile", "unknown")),
        parameter_profile=str(config.get("parameter_profile", "unknown")),
        numerical_profile=str(config.get("numerical_profile", "unknown")),
        has_time_history=bool(histories) and all(path.is_file() for path in histories),
        exact_rng_checkpoint=bool(checkpoints) and all(path.is_file() for path in checkpoints),
        parent_path=(
            continuation.get("parent_run_directory")
            or stack_manifest.get("parent_run_directory")
        ),
    )


def scan_runs(roots: Iterable[str | Path]) -> list[SimulationRecord]:
    """Scan result roots without descending into per-layer directories."""
    records: dict[str, SimulationRecord] = {}
    for root_value in roots:
        root = Path(root_value).expanduser()
        if not root.exists():
            continue
        candidates = [root] if (root / "config.json").is_file() else []
        for stack_file in root.rglob("mechanical_spatial_stack.npz"):
            candidates.append(stack_file.parent)
        for candidate in candidates:
            record = inspect_run(candidate)
            if record is not None:
                records[record.path] = record
    return sorted(records.values(), key=lambda item: item.created_local, reverse=True)
