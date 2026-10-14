"""Extend a completed stack with new independently generated mechanical layers."""

from __future__ import annotations

import json
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from typing import Callable

import numpy as np

from .config import Program19Config
from .mechanical_runner import _write_manifest, run_pipeline
from .persistence import atomic_savez, atomic_write_json
from .reporting import save_stack_slices_figure
from .transport import accessible_centers, classify_void_components


ProgressCallback = Callable[[dict], None]


def referenced_layer_directories(run: str | Path) -> list[Path]:
    run = Path(run).resolve()
    stack_manifest = run / "stack_manifest.json"
    if stack_manifest.is_file():
        data = json.loads(stack_manifest.read_text(encoding="utf-8"))
        result = [Path(value).resolve() for value in data.get("layer_directories", [])]
    else:
        root = run / "mechanical_layers"
        result = sorted(
            path.resolve() for path in root.iterdir()
            if path.is_dir() and path.name.startswith("layer_")
        ) if root.is_dir() else []
    if not result:
        raise RuntimeError(f"No mechanical layers are referenced by {run}")
    return result


def used_layer_seeds(roots: list[str | Path]) -> set[int]:
    seeds: set[int] = set()
    for value in roots:
        root = Path(value).expanduser()
        if not root.exists():
            continue
        for stack_path in root.rglob("mechanical_spatial_stack.npz"):
            try:
                with np.load(stack_path, allow_pickle=False) as saved:
                    if "layer_seeds" in saved.files:
                        seeds.update(int(item) for item in np.asarray(saved["layer_seeds"]).ravel())
            except (OSError, ValueError):
                continue
    return seeds


def _new_directory(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = root / f"program19_3_extended_stack_{timestamp}"
    suffix = 1
    while candidate.exists():
        candidate = root / f"program19_3_extended_stack_{timestamp}_{suffix:02d}"
        suffix += 1
    candidate.mkdir()
    return candidate


def extend_stack(
    source_run: str | Path,
    *,
    target_layer_count: int,
    first_new_seed: int,
    output_root: str | Path,
    mechanical_workers: int | None = None,
    forbidden_seeds: set[int] | None = None,
    progress_callback: ProgressCallback | None = None,
) -> Path:
    """Generate only the added layers and create a non-destructive virtual combined stack."""
    source = Path(source_run).resolve()
    if (source / "continuation_manifest.json").is_file():
        raise ValueError(
            "Extend the original run first, then continue the enlarged stack in time. "
            "Direct extension of a continuation lineage is not yet supported."
        )
    old_layers = referenced_layer_directories(source)
    old_count = len(old_layers)
    if target_layer_count <= old_count:
        raise ValueError(f"target_layer_count must exceed the current {old_count} layers")
    added_count = target_layer_count - old_count
    source_config = Program19Config.load_json(source / "config.json")
    if mechanical_workers is not None:
        source_config = replace(source_config, mechanical_workers=int(mechanical_workers))
        source_config.validate()
    new_seeds = [first_new_seed + index * source_config.layer_seed_stride for index in range(added_count)]
    old_seeds = set()
    for layer in old_layers:
        summary_path = layer / "layer_summary.json"
        if summary_path.is_file():
            old_seeds.add(int(json.loads(summary_path.read_text(encoding="utf-8"))["seed"]))
    collisions = (old_seeds | (forbidden_seeds or set())) & set(new_seeds)
    if collisions:
        raise ValueError(f"New layer seeds collide with existing seeds: {sorted(collisions)}")

    result = _new_directory(Path(output_root).expanduser().resolve())
    component_root = result / "_components"
    component_root.mkdir()
    added_config = replace(
        source_config,
        mechanical_layers=added_count,
        seed=int(first_new_seed),
        output_root=str(component_root),
        run_transport=False,
    )

    def relay(payload: dict) -> None:
        if progress_callback:
            progress_callback({**payload, "stage": "generating_added_layers"})

    _, generated_run, generated_status = run_pipeline(added_config, progress_callback=relay)
    if generated_status != "completed":
        raise RuntimeError(f"Added-layer generation did not complete: {generated_run}")
    new_layers = referenced_layer_directories(generated_run)
    all_layers = old_layers + new_layers
    combined_config = replace(
        source_config, mechanical_layers=target_layer_count, output_root=str(Path(output_root).resolve())
    )
    _write_manifest(result, combined_config)

    solids = []
    regions = []
    layer_seeds = []
    layer_times = []
    extent = None
    dx_um = None
    layer_rows = []
    for index, layer in enumerate(all_layers):
        with np.load(layer / "final_masks.npz", allow_pickle=False) as saved:
            solids.append(np.asarray(saved["solid_mask"], dtype=bool))
            regions.append(np.asarray(saved["analyzed_region_mask"], dtype=bool))
            extent = np.asarray(saved["extent_um"])
            dx_um = float(np.asarray(saved["dx_um"]).item())
            layer_times.append(float(np.asarray(saved["simulated_time_s"]).item()))
        summary = json.loads((layer / "layer_summary.json").read_text(encoding="utf-8"))
        layer_seeds.append(int(summary["seed"]))
        layer_rows.append(
            {
                "layer_index_one_based": index + 1,
                "seed": int(summary["seed"]),
                "simulated_time_s": layer_times[-1],
                "source_layer_directory": str(layer),
                "reused_from_parent": index < old_count,
            }
        )
    stack = np.asarray(solids, dtype=bool)
    region_stack = np.asarray(regions, dtype=bool)
    atomic_savez(
        result / "mechanical_spatial_stack.npz",
        solid_layers=stack,
        analyzed_region_masks=region_stack,
        extent_um=extent,
        dx_um=np.asarray(dx_um),
        dz_um=np.asarray(combined_config.layer_spacing_um),
        layer_z_um=np.arange(target_layer_count, dtype=float) * combined_config.layer_spacing_um,
        layer_seeds=np.asarray(layer_seeds, dtype=int),
        simulated_time_s_by_layer=np.asarray(layer_times),
        source=np.asarray("virtual stack: reused parent layers plus newly generated layers"),
        axis_order=np.asarray("z,y,x"),
    )
    stack_manifest = {
        "program": "Program19.3",
        "operation": "extend_stack_with_independent_mechanical_layers",
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "parent_run_directory": str(source),
        "current_layer_count": old_count,
        "target_layer_count": target_layer_count,
        "added_layer_count": added_count,
        "first_new_seed": first_new_seed,
        "new_seeds": new_seeds,
        "all_seeds": layer_seeds,
        "layer_directories": [str(path) for path in all_layers],
        "generated_component_run": str(generated_run),
        "mechanical_semantics": "all layers are independent 2-D mechanical realizations",
    }
    atomic_write_json(result / "stack_manifest.json", stack_manifest)
    atomic_write_json(result / "layer_index.json", {"layers": layer_rows})

    pore_classes, pore_summary = classify_void_components(stack, region_stack)
    accessible, edt = accessible_centers(
        stack, region_stack, combined_config.particle_radius_um,
        float(dx_um), combined_config.layer_spacing_um,
    )
    atomic_savez(
        result / "pore_analysis.npz",
        pore_class=pore_classes,
        accessible_centers=accessible,
        void_distance_um=edt,
        analyzed_region_masks=region_stack,
        particle_radius_um=np.asarray(combined_config.particle_radius_um),
    )
    pore_summary["accessible_center_voxels"] = int(np.count_nonzero(accessible))
    pore_summary["accessible_porosity_for_particle_centers"] = (
        pore_summary["accessible_center_voxels"] / max(pore_summary["analyzed_region_voxels"], 1)
    )
    pore_summary["mechanical_layers_requested"] = target_layer_count
    pore_summary["mechanical_layers_completed"] = target_layer_count
    atomic_write_json(result / "pore_summary.json", pore_summary)
    atomic_write_json(
        result / "run_metrics.json",
        {
            "mechanical_layers_requested": target_layer_count,
            "mechanical_layers_completed": target_layer_count,
            "geometric_porosity": pore_summary["geometric_porosity"],
            "accessible_porosity_for_particle_centers": pore_summary[
                "accessible_porosity_for_particle_centers"
            ],
            "parent_run_directory": str(source),
            "added_layer_count": added_count,
        },
    )
    save_stack_slices_figure(result / "mechanical_stack_and_pore_slices.png", stack, pore_classes, extent)
    atomic_write_json(
        result / "status.json",
        {
            "status": "completed",
            "requested_mechanical_layers": target_layer_count,
            "completed_or_partial_layers_saved": target_layer_count,
            "error": None,
        },
    )
    if progress_callback:
        progress_callback({"type": "progress", "stage": "complete", "fraction": 1.0})
    return result
