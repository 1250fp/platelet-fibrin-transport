"""Continue completed mechanical layers in time without rebuilding their geometry."""

from __future__ import annotations

import json
import math
import os
import threading
import time
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path
from typing import Callable

import numpy as np

from .config import Program19Config
from .integrator import integrate_interval
from .layers import rasterize_solid
from .mechanical_runner import (
    _filopodia_snapshot,
    _save_intermediate_masks,
    _save_trajectory,
    _write_csv,
    _write_manifest,
)
from .parallel import (
    effective_worker_count,
    numeric_thread_limit,
    run_process_jobs,
    throttled_progress_relay,
)
from .persistence import atomic_savez, atomic_write_json, load_rng, load_state, save_checkpoint
from .reporting import save_mechanics_figure, save_region_figure, save_stack_slices_figure
from .state import SimulationState, StepDiagnostics, state_arrays
from .transport import accessible_centers, classify_void_components


ProgressCallback = Callable[[dict], None]


def _new_directory(root: Path, source: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = f"program19_3_continuation_{source.name[-48:]}_{timestamp}"
    candidate = root / base
    suffix = 1
    while candidate.exists():
        candidate = root / f"{base}_{suffix:02d}"
        suffix += 1
    candidate.mkdir()
    return candidate


def _layer_directories(run: Path) -> list[Path]:
    stack_manifest = run / "stack_manifest.json"
    if stack_manifest.is_file():
        data = json.loads(stack_manifest.read_text(encoding="utf-8"))
        result = [Path(value).resolve() for value in data.get("layer_directories", [])]
        if not result:
            raise RuntimeError(f"No layers are referenced by {stack_manifest}")
        return result
    layers = run / "mechanical_layers"
    if not layers.is_dir():
        raise FileNotFoundError(f"No mechanical_layers folder in {run}")
    result = sorted(path for path in layers.iterdir() if path.is_dir() and path.name.startswith("layer_"))
    if not result:
        raise RuntimeError(f"No saved mechanical layers in {run}")
    return result


def _matching_checkpoint(layer: Path, state_time_s: float) -> Path | None:
    checkpoint = layer / "checkpoint.npz"
    if not checkpoint.is_file():
        return None
    try:
        with np.load(checkpoint, allow_pickle=False) as saved:
            time_s = float(np.asarray(saved["time_s"]).item())
        return checkpoint if math.isclose(time_s, state_time_s, abs_tol=1e-10) else None
    except (OSError, KeyError, ValueError):
        return None


def _continue_one_layer(
    layer_index: int,
    layer_count: int,
    source_layer_text: str,
    config_values: dict,
    start_time: float,
    additional_time_s: float,
    target_time: float,
    layer_dir_text: str,
    progress_callback: ProgressCallback | None,
    pause_event,
    cancel_event,
) -> tuple[SimulationState, str, np.ndarray, np.ndarray, np.ndarray, float, dict]:
    """Continue and persist one saved layer; safe for serial or spawned use."""
    source_layer = Path(source_layer_text)
    layer_dir = Path(layer_dir_text)
    layer_config = Program19Config.from_dict(config_values)
    state = load_state(source_layer / "final_state.npz")
    layer_seed = int(layer_config.seed)
    layer_dir.mkdir(exist_ok=False)
    layer_config.save_json(layer_dir / "config.json")
    atomic_savez(layer_dir / "initial_state.npz", **state_arrays(state))

    fallback = np.random.SeedSequence(
        [
            layer_seed,
            3903,
            int(round(start_time * 1_000_000)),
            int(round(additional_time_s * 1_000_000)),
        ]
    )
    rng, rng_policy = load_rng(
        _matching_checkpoint(source_layer, state.time_s), fallback_seed=fallback
    )
    frame_times = [state.time_s]
    node_frames = [state.nodes_um.copy()]
    platelet_frames = [state.platelets_um.copy()]
    cohesion_frames = [np.asarray(state.cohesion_pairs, dtype=int).reshape((-1, 2))]
    cohesion_reference_frames = [np.asarray(state.cohesion_L0_um, dtype=float)]
    cohesion_birth_frames = [np.asarray(state.cohesion_birth_time_s, dtype=float)]
    filopodia_frames = [_filopodia_snapshot(state)]
    rows: list[dict] = []
    next_output = state.time_s + layer_config.output_interval_s
    layer_wall = time.perf_counter()
    last_diagnostics = StepDiagnostics()
    layer_status = "running"
    layer_error: Exception | None = None
    try:
        while state.time_s < target_time - 1e-15:
            if cancel_event.is_set():
                layer_status = "cancelled"
                break
            while not pause_event.is_set():
                if cancel_event.wait(0.1):
                    layer_status = "cancelled"
                    break
            if layer_status == "cancelled":
                break
            interval = min(layer_config.dt_s, target_time - state.time_s)
            last_diagnostics = integrate_interval(state, layer_config, rng, interval)
            should_store = state.time_s >= next_output - 1e-12 or math.isclose(
                state.time_s, target_time, abs_tol=1e-12
            )
            if should_store:
                rows.append(
                    {
                        "step": state.step,
                        "simulated_time_s": state.time_s,
                        "wall_elapsed_s": time.perf_counter() - layer_wall,
                        **asdict(last_diagnostics),
                    }
                )
                frame_times.append(state.time_s)
                node_frames.append(state.nodes_um.copy())
                platelet_frames.append(state.platelets_um.copy())
                cohesion_frames.append(
                    np.asarray(state.cohesion_pairs, dtype=int).reshape((-1, 2)).copy()
                )
                cohesion_reference_frames.append(
                    np.asarray(state.cohesion_L0_um, dtype=float).copy()
                )
                cohesion_birth_frames.append(
                    np.asarray(state.cohesion_birth_time_s, dtype=float).copy()
                )
                filopodia_frames.append(_filopodia_snapshot(state))
                save_checkpoint(
                    layer_dir / "checkpoint.npz",
                    state,
                    rng,
                    layer_config,
                    source="Program19.3 temporal continuation",
                )
                while next_output <= state.time_s + 1e-12:
                    next_output += layer_config.output_interval_s
            if progress_callback:
                local_fraction = (state.time_s - start_time) / additional_time_s
                progress_callback(
                    {
                        "type": "progress",
                        "stage": "continuing_mechanics",
                        "fraction": (layer_index + local_fraction) / layer_count,
                        "layer_index": layer_index + 1,
                        "layer_count": layer_count,
                        "layer_fraction": local_fraction,
                        "simulated_time_s": state.time_s,
                        "cohesion_links": len(state.cohesion_pairs),
                    }
                )
        if layer_status == "running":
            layer_status = "completed"
    except Exception as exc:
        layer_status = "failed"
        layer_error = exc
    finally:
        atomic_savez(layer_dir / "final_state.npz", **state_arrays(state))
        save_checkpoint(
            layer_dir / "checkpoint.npz",
            state,
            rng,
            layer_config,
            source="Program19.3 temporal continuation final state",
        )
        if layer_config.save_trajectory:
            _save_trajectory(
                layer_dir / "trajectory.npz",
                state,
                layer_config,
                frame_times,
                node_frames,
                platelet_frames,
                cohesion_frames,
                cohesion_reference_frames,
                cohesion_birth_frames,
                filopodia_frames,
            )
        if layer_config.save_intermediate_masks:
            _save_intermediate_masks(
                layer_dir / "intermediate_masks.npz",
                frame_times,
                node_frames,
                platelet_frames,
                state.backbone_pairs,
                layer_config,
            )
        _write_csv(layer_dir / "diagnostics.csv", rows)
        solid, region, extent, dx_um = rasterize_solid(state, layer_config)
        atomic_savez(
            layer_dir / "final_masks.npz",
            solid_mask=solid,
            analyzed_region_mask=region,
            extent_um=np.asarray(extent),
            dx_um=np.asarray(dx_um),
            simulated_time_s=np.asarray(state.time_s),
            source=np.asarray("continued explicit mechanical state"),
        )
        save_mechanics_figure(
            layer_dir / "final_mechanics.png", state, layer_config, layer_status
        )
        save_region_figure(layer_dir / "analyzed_region_2d.png", solid, region, extent)
        summary = {
            "layer_index_one_based": layer_index + 1,
            "seed": layer_seed,
            "status": layer_status,
            "source_layer_directory": str(source_layer),
            "source_time_s": start_time,
            "simulated_time_s": state.time_s,
            "additional_time_completed_s": state.time_s - start_time,
            "completed_outer_steps": state.step,
            "saved_intermediate_states": len(frame_times),
            "node_count": len(state.nodes_um),
            "platelet_count": len(state.platelets_um),
            "cohesion_link_count": len(state.cohesion_pairs),
            "rng_policy": rng_policy,
            "wall_elapsed_s": time.perf_counter() - layer_wall,
            "worker_process_id": os.getpid(),
            "mechanical_workers_requested": layer_config.mechanical_workers,
            "error": str(layer_error) if layer_error else None,
        }
        atomic_write_json(layer_dir / "layer_summary.json", summary)
    if layer_error:
        raise RuntimeError(f"Layer {layer_index + 1} continuation failed: {layer_error}") from layer_error
    return state, layer_status, solid, region, extent, float(dx_um), summary


def _parallel_continuation_worker(
    layer_index: int,
    layer_count: int,
    source_layer_text: str,
    config_values: dict,
    start_time: float,
    additional_time_s: float,
    target_time: float,
    layer_dir_text: str,
    pause_event,
    cancel_event,
    progress_queue,
):
    relay = throttled_progress_relay(progress_queue)
    with numeric_thread_limit():
        return _continue_one_layer(
            layer_index,
            layer_count,
            source_layer_text,
            config_values,
            start_time,
            additional_time_s,
            target_time,
            layer_dir_text,
            relay,
            pause_event,
            cancel_event,
        )


def continue_mechanical_run(
    source_run: str | Path,
    *,
    additional_time_s: float,
    output_root: str | Path,
    mechanical_workers: int | None = None,
    progress_callback: ProgressCallback | None = None,
    pause_event: threading.Event | None = None,
    cancel_event: threading.Event | None = None,
) -> Path:
    """Continue every saved layer by the same additional physical duration."""
    if not np.isfinite(additional_time_s) or additional_time_s <= 0.0:
        raise ValueError("additional_time_s must be positive")
    source = Path(source_run).resolve()
    source_layers = _layer_directories(source)
    source_config = Program19Config.load_json(source / "config.json")
    if mechanical_workers is not None:
        source_config = replace(source_config, mechanical_workers=int(mechanical_workers))
        source_config.validate()
    if len(source_layers) != source_config.mechanical_layers:
        raise RuntimeError(
            f"Expected {source_config.mechanical_layers} complete layers, found {len(source_layers)}"
        )
    start_times = []
    for layer in source_layers:
        with np.load(layer / "final_state.npz", allow_pickle=False) as saved:
            start_times.append(float(np.asarray(saved["time_s"]).item()))
    start_times = np.asarray(start_times, dtype=float)
    if not np.allclose(start_times, start_times[0], atol=1e-10):
        raise RuntimeError("Source layers do not share the same final mechanical time")
    start_time = float(start_times[0])
    target_time = start_time + float(additional_time_s)
    config = replace(source_config, duration_s=target_time, run_transport=False)
    result = _new_directory(Path(output_root).expanduser().resolve(), source)
    _write_manifest(result, config)
    layers_root = result / "mechanical_layers"
    layers_root.mkdir()
    pause_event = pause_event or threading.Event()
    pause_event.set()
    cancel_event = cancel_event or threading.Event()
    continuation_manifest = {
        "program": "Program19.3",
        "operation": "continue_mechanical_evolution",
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "parent_run_directory": str(source),
        "source_time_s": start_time,
        "additional_time_s": float(additional_time_s),
        "target_time_s": target_time,
        "layer_count": len(source_layers),
        "mechanical_workers_requested": config.mechanical_workers,
        "mechanical_workers_effective": effective_worker_count(
            config.mechanical_workers, len(source_layers)
        ),
        "parameter_policy": "mechanical parameters are inherited unchanged; only final time changes",
        "rng_policy_by_layer": [],
    }
    atomic_write_json(result / "continuation_manifest.json", continuation_manifest)

    wall_start = time.perf_counter()
    final_states = []
    solids = []
    regions = []
    summaries = []
    status = "running"
    error: Exception | None = None
    extent = None
    dx_um = None
    try:
        completed_results: dict[int, tuple] = {}

        def accept_result(layer_index: int, item: object) -> None:
            nonlocal extent, dx_um
            completed_results[layer_index] = item  # type: ignore[assignment]
            ordered = [completed_results[index] for index in sorted(completed_results)]
            final_states[:] = [entry[0] for entry in ordered]
            solids[:] = [entry[2] for entry in ordered]
            regions[:] = [entry[3] for entry in ordered]
            summaries[:] = [entry[6] for entry in ordered]
            extent = ordered[-1][4]
            dx_um = ordered[-1][5]
            continuation_manifest["rng_policy_by_layer"] = [
                {
                    "layer": int(summary["layer_index_one_based"]),
                    "seed": int(summary["seed"]),
                    "policy": summary["rng_policy"],
                }
                for summary in summaries
            ]
            atomic_write_json(result / "continuation_manifest.json", continuation_manifest)
            _write_csv(
                result / "layer_index.csv",
                [
                    {
                        **summary,
                        "folder": (
                            f"mechanical_layers/layer_"
                            f"{int(summary['layer_index_one_based']):03d}_seed{summary['seed']}"
                        ),
                    }
                    for summary in summaries
                ],
            )

        jobs = []
        for layer_index, source_layer in enumerate(source_layers):
            source_summary_path = source_layer / "layer_summary.json"
            source_summary = (
                json.loads(source_summary_path.read_text(encoding="utf-8"))
                if source_summary_path.is_file()
                else {}
            )
            layer_seed = int(
                source_summary.get(
                    "seed", source_config.seed + layer_index * source_config.layer_seed_stride
                )
            )
            layer_config = replace(config, seed=layer_seed)
            layer_dir = layers_root / f"layer_{layer_index + 1:03d}_seed{layer_seed}"
            jobs.append(
                (
                    layer_index,
                    len(source_layers),
                    str(source_layer),
                    layer_config.to_dict(),
                    start_time,
                    float(additional_time_s),
                    target_time,
                    str(layer_dir),
                )
            )

        workers = effective_worker_count(config.mechanical_workers, len(source_layers))
        if workers == 1:
            for job in jobs:
                if cancel_event.is_set():
                    status = "cancelled"
                    break
                item = _continue_one_layer(
                    *job,
                    progress_callback,
                    pause_event,
                    cancel_event,
                )
                accept_result(int(job[0]), item)
                if item[1] == "cancelled":
                    status = "cancelled"
                    break
        else:
            run_process_jobs(
                _parallel_continuation_worker,
                jobs,
                requested_workers=config.mechanical_workers,
                layer_count=len(source_layers),
                stage="continuing_mechanics",
                progress_callback=progress_callback,
                result_callback=accept_result,
                pause_event=pause_event,
                cancel_event=cancel_event,
            )
            if cancel_event.is_set() or any(
                summary.get("status") == "cancelled" for summary in summaries
            ) or len(summaries) < len(source_layers):
                status = "cancelled"

        if status == "running":
            status = "completed"
        if not solids:
            raise RuntimeError("No continued layer reached a saveable state")
        stack = np.asarray(solids, dtype=bool)
        region_stack = np.asarray(regions, dtype=bool)
        layer_seeds = np.asarray([item["seed"] for item in summaries], dtype=int)
        layer_times = np.asarray([item["simulated_time_s"] for item in summaries], dtype=float)
        atomic_savez(
            result / "mechanical_spatial_stack.npz",
            solid_layers=stack,
            analyzed_region_masks=region_stack,
            extent_um=np.asarray(extent),
            dx_um=np.asarray(dx_um),
            dz_um=np.asarray(config.layer_spacing_um),
            layer_z_um=np.arange(len(stack), dtype=float) * config.layer_spacing_um,
            layer_seeds=layer_seeds,
            simulated_time_s_by_layer=layer_times,
            source=np.asarray("continued independent explicit mechanical layers"),
            axis_order=np.asarray("z,y,x"),
        )
        pore_classes, pore_summary = classify_void_components(stack, region_stack)
        accessible, edt = accessible_centers(
            stack, region_stack, config.particle_radius_um, float(dx_um), config.layer_spacing_um
        )
        atomic_savez(
            result / "pore_analysis.npz",
            pore_class=pore_classes,
            accessible_centers=accessible,
            void_distance_um=edt,
            analyzed_region_masks=region_stack,
            particle_radius_um=np.asarray(config.particle_radius_um),
        )
        pore_summary["accessible_center_voxels"] = int(np.count_nonzero(accessible))
        pore_summary["accessible_porosity_for_particle_centers"] = (
            pore_summary["accessible_center_voxels"]
            / max(pore_summary["analyzed_region_voxels"], 1)
        )
        pore_summary["mechanical_layers_requested"] = config.mechanical_layers
        pore_summary["mechanical_layers_completed"] = len(stack)
        atomic_write_json(result / "pore_summary.json", pore_summary)
        save_stack_slices_figure(result / "mechanical_stack_and_pore_slices.png", stack, pore_classes, extent)
        atomic_write_json(
            result / "run_metrics.json",
            {
                "mechanical_layers_requested": config.mechanical_layers,
                "mechanical_layers_completed": len(stack),
                "mechanical_workers_requested": config.mechanical_workers,
                "mechanical_workers_effective": effective_worker_count(
                    config.mechanical_workers, len(source_layers)
                ),
                "source_time_s": start_time,
                "final_time_s": float(np.min(layer_times)),
                "additional_time_requested_s": additional_time_s,
                "geometric_porosity": pore_summary["geometric_porosity"],
                "accessible_porosity_for_particle_centers": pore_summary[
                    "accessible_porosity_for_particle_centers"
                ],
                "autoregression_used": False,
                "interlayer_force_coupling": False,
            },
        )
        if progress_callback:
            progress_callback({"type": "progress", "stage": "complete", "fraction": 1.0})
    except Exception as exc:
        status = "failed"
        error = exc
    finally:
        atomic_write_json(
            result / "status.json",
            {
                "status": status,
                "requested_mechanical_layers": config.mechanical_layers,
                "completed_or_partial_layers_saved": len(solids),
                "mechanical_workers_requested": config.mechanical_workers,
                "mechanical_workers_effective": effective_worker_count(
                    config.mechanical_workers, len(source_layers)
                ),
                "wall_elapsed_s": time.perf_counter() - wall_start,
                "error": str(error) if error else None,
            },
        )
    if error:
        raise RuntimeError(f"{error} Partial continuation saved in {result}") from error
    return result
