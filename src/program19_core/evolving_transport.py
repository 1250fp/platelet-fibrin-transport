"""Central-release first passage through a mechanically evolving saved stack."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable

import numpy as np

from .persistence import atomic_savez, atomic_write_json
from .transport import accessible_centers
from .transport_steps import (
    ACCEPT, ESCAPE_Z_MIN, ESCAPE_Z_MAX, REJECT_SOLID, REJECT_LATERAL,
    TRANSPORT_KERNEL_VERSION, COLLISION_RULE, bounded_step_duration,
    empirical_survival, trace_proposals,
)

try:
    from scipy.ndimage import distance_transform_edt
except ImportError as exc:  # pragma: no cover
    distance_transform_edt = None
    SCIPY_ERROR = exc
else:
    SCIPY_ERROR = None


ProgressCallback = Callable[[dict], None]


@dataclass(frozen=True)
class StackHistory:
    times_s: np.ndarray
    solid_masks: np.ndarray  # time,z,y,x
    analyzed_regions: np.ndarray  # time,z,y,x
    dx_um: float
    dz_um: float
    lineage: tuple[str, ...]


def _continuation_parent(run: Path) -> Path | None:
    manifest_path = run / "continuation_manifest.json"
    if not manifest_path.is_file():
        return None
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    value = data.get("parent_run_directory")
    return Path(value).resolve() if value else None


def _run_lineage(run: Path, seen: set[Path] | None = None) -> list[Path]:
    seen = seen or set()
    resolved = run.resolve()
    if resolved in seen:
        raise RuntimeError(f"Continuation lineage cycle at {resolved}")
    seen.add(resolved)
    parent = _continuation_parent(resolved)
    return (_run_lineage(parent, seen) if parent is not None else []) + [resolved]


def _layer_directories(run: Path) -> list[Path]:
    stack_manifest = run / "stack_manifest.json"
    if stack_manifest.is_file():
        data = json.loads(stack_manifest.read_text(encoding="utf-8"))
        layers = [Path(value).resolve() for value in data.get("layer_directories", [])]
        if not layers:
            raise RuntimeError(f"No layers are referenced by {stack_manifest}")
        return layers
    root = run / "mechanical_layers"
    if not root.is_dir():
        raise FileNotFoundError(f"Mechanical-layer directory not found: {root}")
    layers = sorted(path for path in root.iterdir() if path.is_dir() and path.name.startswith("layer_"))
    if not layers:
        raise RuntimeError(f"No mechanical layers found in {run}")
    return layers


def _load_segment(run: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    config = json.loads((run / "config.json").read_text(encoding="utf-8"))
    dz_um = float(config["layer_spacing_um"])
    layer_times: list[np.ndarray] = []
    layer_solids: list[np.ndarray] = []
    layer_regions: list[np.ndarray] = []
    dx_values: list[float] = []
    for layer in _layer_directories(run):
        history_path = layer / "intermediate_masks.npz"
        if history_path.is_file():
            with np.load(history_path, allow_pickle=False) as saved:
                layer_times.append(np.asarray(saved["times_s"], dtype=float))
                layer_solids.append(np.asarray(saved["solid_masks"], dtype=bool))
                layer_regions.append(np.asarray(saved["analyzed_region_masks"], dtype=bool))
                dx_values.append(float(np.asarray(saved["dx_um"]).item()))
        else:
            final_path = layer / "final_masks.npz"
            if not final_path.is_file():
                raise FileNotFoundError(f"No mask history or final mask in {layer}")
            with np.load(final_path, allow_pickle=False) as saved:
                layer_times.append(np.asarray([float(saved["simulated_time_s"].item())]))
                layer_solids.append(np.asarray(saved["solid_mask"], dtype=bool)[None, ...])
                layer_regions.append(
                    np.asarray(saved["analyzed_region_mask"], dtype=bool)[None, ...]
                )
                dx_values.append(float(np.asarray(saved["dx_um"]).item()))

    reference_times = layer_times[0]
    reference_shape = layer_solids[0].shape
    for index, (times, solids, regions) in enumerate(
        zip(layer_times, layer_solids, layer_regions), start=1
    ):
        if solids.shape != reference_shape or regions.shape != reference_shape:
            raise ValueError(f"Layer {index} has a different time-mask shape")
        if times.shape != reference_times.shape or not np.allclose(times, reference_times, atol=1e-10):
            raise ValueError(f"Layer {index} does not share the common mechanical time grid")
    if not np.allclose(dx_values, dx_values[0]):
        raise ValueError("Layer raster spacings differ")
    # Older Program19.2 runs can contain the numerically reached end time and
    # the exact requested end time separated only by round-off. Keep the last
    # of such adjacent frames so the mechanical clock remains strictly increasing.
    keep = np.r_[np.diff(reference_times) > 1e-10, True]
    solids = np.stack([item[keep] for item in layer_solids], axis=1)
    regions = np.stack([item[keep] for item in layer_regions], axis=1)
    return reference_times[keep], solids, regions, dx_values[0], dz_um


def load_stack_history(run_directory: str | Path, include_lineage: bool = True) -> StackHistory:
    run = Path(run_directory).resolve()
    chain = _run_lineage(run) if include_lineage else [run]
    times_out: list[np.ndarray] = []
    solids_out: list[np.ndarray] = []
    regions_out: list[np.ndarray] = []
    dx_um = dz_um = None
    last_time = -math.inf
    expected_spatial_shape = None
    for segment in chain:
        times, solids, regions, current_dx, current_dz = _load_segment(segment)
        if dx_um is None:
            dx_um, dz_um = current_dx, current_dz
            expected_spatial_shape = solids.shape[1:]
        elif (
            not math.isclose(dx_um, current_dx)
            or not math.isclose(dz_um, current_dz)
            or solids.shape[1:] != expected_spatial_shape
        ):
            raise ValueError("Continuation segments do not share the same stack discretization")
        keep = times > last_time + 1e-10
        if np.any(keep):
            times_out.append(times[keep])
            solids_out.append(solids[keep])
            regions_out.append(regions[keep])
            last_time = float(times[keep][-1])
    if not times_out:
        raise RuntimeError("No time-resolved stack masks were found")
    return StackHistory(
        times_s=np.concatenate(times_out),
        solid_masks=np.concatenate(solids_out, axis=0),
        analyzed_regions=np.concatenate(regions_out, axis=0),
        dx_um=float(dx_um),
        dz_um=float(dz_um),
        lineage=tuple(str(path) for path in chain),
    )


class _AccessibilityCache:
    def __init__(self, history: StackHistory, radius_um: float):
        self.history = history
        self.radius_um = float(radius_um)
        self.values: dict[int, np.ndarray] = {}

    def get(self, frame: int) -> np.ndarray:
        frame = int(frame)
        if frame not in self.values:
            accessible, _ = accessible_centers(
                self.history.solid_masks[frame],
                self.history.analyzed_regions[frame],
                self.radius_um,
                self.history.dx_um,
                self.history.dz_um,
            )
            self.values[frame] = accessible
        return self.values[frame]


def _frame_at(times: np.ndarray, time_s: float) -> int:
    return int(np.clip(np.searchsorted(times, time_s, side="right") - 1, 0, len(times) - 1))


def _voxel_indices(positions: np.ndarray, spacing: np.ndarray) -> np.ndarray:
    return np.floor(positions / spacing).astype(int)


def _is_accessible(accessible: np.ndarray, positions: np.ndarray, spacing: np.ndarray) -> np.ndarray:
    indices = _voxel_indices(positions, spacing)
    shape = np.asarray(accessible.shape)
    inside = np.all((indices >= 0) & (indices < shape), axis=1)
    result = np.zeros(len(positions), dtype=bool)
    valid = np.flatnonzero(inside)
    if len(valid):
        idx = indices[valid]
        result[valid] = accessible[idx[:, 0], idx[:, 1], idx[:, 2]]
    return result


def _relocate_from_moving_solid(
    positions: np.ndarray,
    active_ids: np.ndarray,
    accessible: np.ndarray,
    spacing: np.ndarray,
) -> int:
    """Project centers newly covered by a moving mask to the nearest accessible voxel."""
    allowed = _is_accessible(accessible, positions[active_ids], spacing)
    invalid_ids = active_ids[~allowed]
    if not len(invalid_ids):
        return 0
    if distance_transform_edt is None:
        raise ImportError(f"SciPy is required for moving-boundary relocation: {SCIPY_ERROR}")
    if not np.any(accessible):
        raise RuntimeError("The evolving stack has no accessible particle-center voxel")
    _, nearest = distance_transform_edt(
        ~accessible, sampling=tuple(spacing), return_indices=True
    )
    indices = _voxel_indices(positions[invalid_ids], spacing)
    indices = np.clip(indices, 0, np.asarray(accessible.shape) - 1)
    target = nearest[:, indices[:, 0], indices[:, 1], indices[:, 2]].T
    positions[invalid_ids] = (target.astype(float) + 0.5) * spacing
    return int(len(invalid_ids))


def time_dependent_first_passage(
    history: StackHistory,
    *,
    release_time_s: float,
    observation_time_s: float,
    particle_radius_um: float,
    diffusion_um2_s: float,
    particle_count: int,
    dt_s: float,
    rng: np.random.Generator,
    geometry_mode: str = "evolving",
    recorded_particles: int = 25,
    trajectory_interval_s: float = 0.01,
    progress_callback: ProgressCallback | None = None,
) -> tuple[dict, dict[str, np.ndarray]]:
    """Run central-release FPT with frozen or piecewise-constant evolving masks."""
    if geometry_mode not in {"evolving", "frozen"}:
        raise ValueError("geometry_mode must be 'evolving' or 'frozen'")
    for name, value in (
        ("observation_time_s", observation_time_s),
        ("particle_radius_um", particle_radius_um),
        ("diffusion_um2_s", diffusion_um2_s),
        ("dt_s", dt_s),
        ("trajectory_interval_s", trajectory_interval_s),
    ):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive")
    if particle_count < 1:
        raise ValueError("particle_count must be at least 1")
    if release_time_s < history.times_s[0] - 1e-10:
        raise ValueError("release_time_s precedes the saved mechanical history")
    if release_time_s > history.times_s[-1] + 1e-10:
        raise ValueError("release_time_s follows the saved mechanical history")
    if geometry_mode == "evolving" and release_time_s + observation_time_s > history.times_s[-1] + 1e-10:
        raise ValueError(
            "The evolving mechanical history ends before release_time_s + observation_time_s"
        )

    cache = _AccessibilityCache(history, particle_radius_um)
    release_frame = _frame_at(history.times_s, release_time_s)
    current_frame = release_frame
    accessible = cache.get(release_frame)
    nz, ny, nx = accessible.shape
    zz, yy, xx = np.indices(accessible.shape)
    central = (
        accessible
        & (zz >= int(0.4 * nz))
        & (zz < max(int(0.6 * nz), int(0.4 * nz) + 1))
        & ((xx - 0.5 * (nx - 1)) ** 2 + (yy - 0.5 * (ny - 1)) ** 2
           <= (0.20 * min(nx, ny)) ** 2)
    )
    starts = np.argwhere(central)
    initialization = "central accessible region"
    if not len(starts):
        starts = np.argwhere(accessible)
        initialization = "entire accessible domain (central region was empty)"
    if not len(starts):
        raise RuntimeError("No accessible voxel is available at the release time")

    selected = starts[rng.integers(0, len(starts), size=particle_count)]
    spacing = np.asarray((history.dz_um, history.dx_um, history.dx_um), dtype=float)
    positions = (selected.astype(float) + 0.5) * spacing
    active = np.ones(particle_count, dtype=bool)
    outcomes = np.zeros(particle_count, dtype=np.uint8)  # 0 censored, 1 z-min, 2 z-max
    escape_times = np.full(particle_count, np.nan)
    step_sigma = math.sqrt(2.0 * diffusion_um2_s * dt_s)
    total_steps = int(math.ceil(observation_time_s / dt_s))
    moving_relocations = 0
    rejected_solid = 0
    rejected_lateral = 0

    recorded = min(max(int(recorded_particles), 0), particle_count)
    record_every = max(1, int(round(trajectory_interval_s / dt_s)))
    trajectory_times = [0.0]
    trajectory_positions = [positions[:recorded].copy()]
    trajectory_outcomes = [outcomes[:recorded].copy()]

    for step in range(1, total_steps + 1):
        h = bounded_step_duration(step, dt_s, observation_time_s)
        if h <= 0:
            break
        relative_time = min(step * dt_s, observation_time_s)
        absolute_time = release_time_s + (step - 1) * dt_s
        if geometry_mode == "evolving":
            frame = _frame_at(history.times_s, absolute_time)
            if frame != current_frame:
                current_frame = frame
                accessible = cache.get(frame)
                moving_relocations += _relocate_from_moving_solid(
                    positions, np.flatnonzero(active), accessible, spacing
                )

        ids = np.flatnonzero(active)
        if not len(ids):
            break
        # Draw for every particle so frozen/evolving paired runs share each
        # particle's Brownian increments even after their active sets diverge.
        noise = rng.normal(0.0, math.sqrt(2.0*diffusion_um2_s*h), size=(particle_count, 3))
        proposal = positions[ids] + noise[ids]
        events, fractions = trace_proposals(accessible, positions[ids], proposal, spacing)
        lower_escape = events == ESCAPE_Z_MIN
        upper_escape = events == ESCAPE_Z_MAX
        escaped = lower_escape | upper_escape
        if np.any(escaped):
            escaped_ids = ids[escaped]
            escape_times[escaped_ids] = relative_time
            outcomes[ids[lower_escape]] = 1
            outcomes[ids[upper_escape]] = 2
            # Record the geometric exit location on the face, not an arbitrary
            # overshoot. The event time still has the transport-step resolution.
            positions[escaped_ids] += (
                proposal[escaped]-positions[escaped_ids]
            )*fractions[escaped,None]
            active[escaped_ids] = False
        allowed = events == ACCEPT
        positions[ids[allowed]] = proposal[allowed]
        rejected_lateral += int(np.count_nonzero(events == REJECT_LATERAL))
        rejected_solid += int(np.count_nonzero(events == REJECT_SOLID))

        if step % record_every == 0 or step == total_steps or not np.any(active):
            trajectory_times.append(relative_time)
            trajectory_positions.append(positions[:recorded].copy())
            trajectory_outcomes.append(outcomes[:recorded].copy())
        if progress_callback and (step % max(1, total_steps // 100) == 0 or step == total_steps):
            progress_callback(
                {
                    "type": "progress",
                    "stage": f"{geometry_mode}_transport",
                    "fraction": step / total_steps,
                    "simulated_time_s": relative_time,
                    "active_particles": int(np.count_nonzero(active)),
                }
            )

    escaped = outcomes > 0
    observed = np.where(escaped, escape_times, observation_time_s)
    survival_grid = np.linspace(0.0, observation_time_s, 101)
    survival = empirical_survival(escape_times, survival_grid, observation_time_s)
    summary = {
        "geometry_mode": geometry_mode,
        "particle_count": int(particle_count),
        "escaped_count": int(np.count_nonzero(escaped)),
        "escaped_z_min_count": int(np.count_nonzero(outcomes == 1)),
        "escaped_z_max_count": int(np.count_nonzero(outcomes == 2)),
        "censored_count": int(np.count_nonzero(~escaped)),
        "escape_fraction": float(np.mean(escaped)),
        "retained_fraction_at_censor": float(np.mean(~escaped)),
        "restricted_mean_first_passage_time_s": float(np.mean(observed)),
        "conditional_mean_escape_time_s": (
            float(np.mean(escape_times[escaped])) if np.any(escaped) else None
        ),
        "release_time_s": float(release_time_s),
        "observation_time_s": float(observation_time_s),
        "mechanical_frame_at_release": int(release_frame),
        "mechanical_frame_time_s": float(history.times_s[release_frame]),
        "mechanical_frame_count": int(len(history.times_s)),
        "mechanical_frame_interval_interpretation": "piecewise-constant saved masks",
        "initialization": initialization,
        "particle_radius_um": float(particle_radius_um),
        "diffusion_um2_s": float(diffusion_um2_s),
        "transport_dt_s": float(dt_s),
        "step_sigma_per_axis_um": float(step_sigma),
        "moving_boundary_relocations": int(moving_relocations),
        "solid_collision_rejections": int(rejected_solid),
        "lateral_boundary_rejections": int(rejected_lateral),
        "transport_kernel_version": TRANSPORT_KERNEL_VERSION,
        "collision_rule": COLLISION_RULE,
        "last_step_rule": "variance 2*D*min(dt, remaining observation time)",
        "exit_definition": "unobstructed first particle-center contact with either z face",
        "escape_time_resolution": "end of transport step; not exact Brownian hitting time",
        "survival_definition": "mean(censored | (observed escape time > t)), 0<=t<=Tc",
        "moving_boundary_rule": (
            "a center newly covered by a changed mask is projected to the nearest "
            "accessible-center voxel; this is steric displacement, not biochemical capture"
        ),
        "outcome_codes": {"0": "retained_at_censor", "1": "escaped_z_min", "2": "escaped_z_max"},
    }
    arrays = {
        "escape_times_s": escape_times,
        "outcomes": outcomes,
        "survival_time_s": survival_grid,
        "survival_probability": survival,
        "trajectory_time_since_release_s": np.asarray(trajectory_times),
        "trajectory_positions_zyx_um": np.asarray(trajectory_positions),
        "trajectory_outcomes": np.asarray(trajectory_outcomes),
        "recorded_particle_ids": np.arange(recorded, dtype=int),
    }
    return summary, arrays


def run_frozen_transport_analysis(
    source_run: str | Path,
    *,
    output_root: str | Path,
    release_time_s: float,
    observation_time_s: float,
    particle_radius_um: float,
    diffusion_um2_s: float,
    particle_count: int,
    dt_s: float,
    seed: int,
    recorded_particles: int = 25,
    trajectory_interval_s: float = 0.01,
    progress_callback: ProgressCallback | None = None,
) -> Path:
    """Run post-hoc transport on one frozen saved mechanical frame.

    The transport observation horizon is independent of the duration of the
    saved mechanical history because the selected geometry remains fixed.
    """
    source = Path(source_run).resolve()
    history = load_stack_history(source, include_lineage=True)
    release_frame = _frame_at(history.times_s, release_time_s)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    source_tag = source.name[-48:]
    result = Path(output_root).expanduser().resolve() / (
        f"program19_3_frozen_transport_{source_tag}_{timestamp}"
    )
    result.mkdir(parents=True, exist_ok=False)
    request = {
        "program": "Program19.3",
        "analysis": "post-hoc central-release transport through a frozen saved mechanical stack",
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "source_run_directory": str(source),
        "source_lineage": list(history.lineage),
        "requested_release_time_s": float(release_time_s),
        "selected_mechanical_frame_index": int(release_frame),
        "selected_mechanical_frame_time_s": float(history.times_s[release_frame]),
        "observation_time_s": float(observation_time_s),
        "observation_horizon_interpretation": (
            "independent of the mechanical-history end because the selected stack is frozen"
        ),
        "particle_radius_um": float(particle_radius_um),
        "diffusion_um2_s": float(diffusion_um2_s),
        "particle_count": int(particle_count),
        "transport_dt_s": float(dt_s),
        "transport_seed": int(seed),
        "recorded_particles": int(recorded_particles),
        "trajectory_interval_s": float(trajectory_interval_s),
        "mechanical_history_start_s": float(history.times_s[0]),
        "mechanical_history_end_s": float(history.times_s[-1]),
        "mechanical_frames": int(len(history.times_s)),
        "stack_layers": int(history.solid_masks.shape[1]),
        "dx_um": history.dx_um,
        "dz_um": history.dz_um,
    }
    atomic_write_json(result / "analysis_manifest.json", request)

    rng = np.random.default_rng(np.random.SeedSequence([int(seed), 9103]))
    summary, arrays = time_dependent_first_passage(
        history,
        release_time_s=release_time_s,
        observation_time_s=observation_time_s,
        particle_radius_um=particle_radius_um,
        diffusion_um2_s=diffusion_um2_s,
        particle_count=particle_count,
        dt_s=dt_s,
        rng=rng,
        geometry_mode="frozen",
        recorded_particles=recorded_particles,
        trajectory_interval_s=trajectory_interval_s,
        progress_callback=progress_callback,
    )
    atomic_savez(result / "frozen_first_passage.npz", **arrays)
    atomic_write_json(result / "frozen_summary.json", summary)

    try:
        import matplotlib.pyplot as plt

        fig, axis = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)
        axis.plot(
            arrays["survival_time_s"], arrays["survival_probability"],
            color="#335c81", linewidth=2,
        )
        axis.set(
            xlabel="Time since release (s)",
            ylabel="Survival probability",
            ylim=(-0.02, 1.02),
            title=f"Frozen mechanical frame at t = {history.times_s[release_frame]:g} s",
        )
        axis.grid(alpha=0.25)
        fig.savefig(result / "survival_frozen.png", dpi=220)
        plt.close(fig)
    except Exception as exc:  # plotting must not invalidate numerical output
        atomic_write_json(result / "plot_warning.json", {"warning": str(exc)})
    return result


def run_evolving_transport_comparison(
    source_run: str | Path,
    *,
    output_root: str | Path,
    release_time_s: float,
    observation_time_s: float,
    particle_radius_um: float,
    diffusion_um2_s: float,
    particle_count: int,
    dt_s: float,
    seed: int,
    recorded_particles: int = 25,
    trajectory_interval_s: float = 0.01,
    include_frozen: bool = True,
    progress_callback: ProgressCallback | None = None,
) -> Path:
    source = Path(source_run).resolve()
    history = load_stack_history(source, include_lineage=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    source_tag = source.name[-48:]
    result = Path(output_root).expanduser().resolve() / (
        f"program19_3_transport_{source_tag}_{timestamp}"
    )
    result.mkdir(parents=True, exist_ok=False)
    request = {
        "program": "Program19.3",
        "analysis": "central-release transport through a time-dependent mechanical stack",
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "source_run_directory": str(source),
        "source_lineage": list(history.lineage),
        "release_time_s": release_time_s,
        "observation_time_s": observation_time_s,
        "particle_radius_um": particle_radius_um,
        "diffusion_um2_s": diffusion_um2_s,
        "particle_count": particle_count,
        "transport_dt_s": dt_s,
        "transport_seed": seed,
        "recorded_particles": recorded_particles,
        "trajectory_interval_s": trajectory_interval_s,
        "include_frozen_comparison": include_frozen,
        "mechanical_history_start_s": float(history.times_s[0]),
        "mechanical_history_end_s": float(history.times_s[-1]),
        "mechanical_frames": int(len(history.times_s)),
        "stack_layers": int(history.solid_masks.shape[1]),
        "dx_um": history.dx_um,
        "dz_um": history.dz_um,
    }
    atomic_write_json(result / "analysis_manifest.json", request)

    modes = ["evolving"] + (["frozen"] if include_frozen else [])
    summaries: dict[str, dict] = {}
    for mode_index, mode in enumerate(modes):
        rng = np.random.default_rng(np.random.SeedSequence([int(seed), 9103]))

        def relay(payload: dict) -> None:
            if progress_callback:
                progress_callback(
                    {
                        **payload,
                        "fraction": (mode_index + float(payload.get("fraction", 0.0))) / len(modes),
                    }
                )

        summary, arrays = time_dependent_first_passage(
            history,
            release_time_s=release_time_s,
            observation_time_s=observation_time_s,
            particle_radius_um=particle_radius_um,
            diffusion_um2_s=diffusion_um2_s,
            particle_count=particle_count,
            dt_s=dt_s,
            rng=rng,
            geometry_mode=mode,
            recorded_particles=recorded_particles,
            trajectory_interval_s=trajectory_interval_s,
            progress_callback=relay,
        )
        summaries[mode] = summary
        atomic_savez(result / f"{mode}_first_passage.npz", **arrays)
        atomic_write_json(result / f"{mode}_summary.json", summary)

    comparison = {
        "definitions": {
            "frozen": "the release-time stack remains fixed during every tracer trajectory",
            "evolving": "saved mechanical masks change with absolute mechanical time",
            "restricted_mean_first_passage_time": "mean(min(T_escape,T_censor))",
        },
        "summaries": summaries,
    }
    if "frozen" in summaries:
        comparison["paired_differences_evolving_minus_frozen"] = {
            "escape_fraction": summaries["evolving"]["escape_fraction"]
            - summaries["frozen"]["escape_fraction"],
            "restricted_mean_first_passage_time_s": summaries["evolving"][
                "restricted_mean_first_passage_time_s"
            ]
            - summaries["frozen"]["restricted_mean_first_passage_time_s"],
        }
    atomic_write_json(result / "comparison_summary.json", comparison)

    try:
        import matplotlib.pyplot as plt

        fig, axis = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)
        for mode in modes:
            with np.load(result / f"{mode}_first_passage.npz", allow_pickle=False) as saved:
                axis.plot(
                    saved["survival_time_s"], saved["survival_probability"],
                    label=mode.capitalize(), linewidth=2,
                )
        axis.set(xlabel="Time since release (s)", ylabel="Survival probability", ylim=(-0.02, 1.02))
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
        fig.savefig(result / "survival_comparison.png", dpi=220)
        plt.close(fig)
    except Exception as exc:  # plotting must not invalidate numerical output
        atomic_write_json(result / "plot_warning.json", {"warning": str(exc)})
    return result
