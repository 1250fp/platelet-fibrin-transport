"""3-D void classification, accessibility, and censored first-passage analysis."""

from __future__ import annotations

import numpy as np

from .transport_steps import (
    ACCEPT, ESCAPE_Z_MIN, ESCAPE_Z_MAX, REJECT_SOLID, REJECT_LATERAL,
    TRANSPORT_KERNEL_VERSION, COLLISION_RULE, bounded_step_duration,
    empirical_survival, trace_proposals,
)

try:
    from scipy.ndimage import distance_transform_edt, generate_binary_structure, label
except ImportError as exc:  # pragma: no cover
    distance_transform_edt = generate_binary_structure = label = None
    SCIPY_ERROR = exc
else:
    SCIPY_ERROR = None


PORE_OUTSIDE_OR_SOLID = 0
PORE_THROUGH_STACK = 1
PORE_BLIND_Z_CONNECTED = 2
PORE_NO_Z_FACE_CONTACT = 3


def _require_scipy() -> None:
    if label is None:
        raise ImportError(f"SciPy is required for pore analysis: {SCIPY_ERROR}")


def _region_stack(analyzed_region: np.ndarray, solid_shape: tuple[int, ...]) -> np.ndarray:
    region = np.asarray(analyzed_region, dtype=bool)
    if region.ndim == 2:
        return np.broadcast_to(region, solid_shape)
    if region.shape != solid_shape:
        raise ValueError(
            f"Analyzed region shape {region.shape} does not match solid stack {solid_shape}"
        )
    return region


def classify_void_components(
    solid_stack: np.ndarray, analyzed_region: np.ndarray
) -> tuple[np.ndarray, dict]:
    """Classify connected void components by contact with the two z faces.

    Class 3 intentionally means *no z-face contact*, not a claim that a cavity
    is completely enclosed with respect to every possible boundary.
    """
    _require_scipy()
    solid = np.asarray(solid_stack, dtype=bool)
    region = _region_stack(analyzed_region, solid.shape)
    void = region & ~solid
    components, count = label(void, structure=generate_binary_structure(3, 1))
    classified = np.zeros_like(components, dtype=np.uint8)
    counts = {
        "through_stack_void_voxels": 0,
        "blind_z_connected_void_voxels": 0,
        "internal_no_z_face_contact_voxels": 0,
        "void_voxels": int(np.count_nonzero(void)),
        "solid_voxels": int(np.count_nonzero(solid & region)),
        "analyzed_region_voxels": int(np.count_nonzero(region)),
        "component_count": int(count),
    }
    for component_id in range(1, count + 1):
        component = components == component_id
        touches_first = bool(np.any(component[0]))
        touches_last = bool(np.any(component[-1]))
        voxel_count = int(np.count_nonzero(component))
        if touches_first and touches_last:
            code = PORE_THROUGH_STACK
            counts["through_stack_void_voxels"] += voxel_count
        elif touches_first or touches_last:
            code = PORE_BLIND_Z_CONNECTED
            counts["blind_z_connected_void_voxels"] += voxel_count
        else:
            code = PORE_NO_Z_FACE_CONTACT
            counts["internal_no_z_face_contact_voxels"] += voxel_count
        classified[component] = code
    total_void = max(counts["void_voxels"], 1)
    total_region = max(counts["analyzed_region_voxels"], 1)
    counts["geometric_porosity"] = counts["void_voxels"] / total_region
    counts["solid_fraction"] = counts["solid_voxels"] / total_region
    counts["through_stack_fraction_of_void"] = counts["through_stack_void_voxels"] / total_void
    counts["blind_z_connected_fraction_of_void"] = counts["blind_z_connected_void_voxels"] / total_void
    counts["internal_no_z_face_contact_fraction_of_void"] = (
        counts["internal_no_z_face_contact_voxels"] / total_void
    )
    return classified, counts


def accessible_centers(
    solid_stack: np.ndarray,
    analyzed_region: np.ndarray,
    particle_radius_um: float,
    dx_um: float,
    dz_um: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Minkowski erosion of void space through an anisotropic physical EDT."""
    _require_scipy()
    solid = np.asarray(solid_stack, dtype=bool)
    region = _region_stack(analyzed_region, solid.shape)
    void = region & ~solid
    distance = distance_transform_edt(void, sampling=(dz_um, dx_um, dx_um))
    return void & (distance >= particle_radius_um), distance


def first_passage_simulation(
    accessible: np.ndarray,
    dx_um: float,
    dz_um: float,
    diffusion_um2_s: float,
    particle_count: int,
    dt_s: float,
    max_time_s: float,
    rng: np.random.Generator,
    injection_rate_cells_s: float = 0.0,
    *,
    return_outcomes: bool = False,
) -> tuple:
    """Off-lattice Brownian first passage through either z face.

    Segment-supercover collisions reject the entire proposal and consume time.
    An unobstructed first z-face contact is an escape, recorded at step end.
    Particles that remain at max_time are right-censored. By default the return
    remains (summary, escape_times, survival_curve); return_outcomes=True adds
    a fourth array (0=censored, 1=z-min, 2=z-max) for shared analysis helpers.
    """
    accessible = np.asarray(accessible, dtype=bool)
    if accessible.ndim != 3 or any(s < 1 for s in accessible.shape):
        raise ValueError("accessible must be a nonempty 3-D mask")
    for name, value in (("dx_um", dx_um), ("dz_um", dz_um), ("dt_s", dt_s),
                        ("max_time_s", max_time_s)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    if not np.isfinite(diffusion_um2_s) or diffusion_um2_s < 0:
        raise ValueError("diffusion_um2_s must be nonnegative and finite")
    if particle_count < 1 or int(particle_count) != particle_count:
        raise ValueError("particle_count must be a positive integer")
    nz, ny, nx = accessible.shape
    zz, yy, xx = np.indices(accessible.shape)
    central = (
        accessible
        & (zz >= int(0.4 * nz))
        & (zz < max(int(0.6 * nz), int(0.4 * nz) + 1))
        & ((xx - 0.5 * (nx - 1)) ** 2 + (yy - 0.5 * (ny - 1)) ** 2 <= (0.20 * min(nx, ny)) ** 2)
    )
    starts = np.argwhere(central)
    initialization = "central accessible region"
    if not len(starts):
        starts = np.argwhere(accessible)
        initialization = "all accessible voxels (central region empty)"
    if not len(starts):
        raise RuntimeError("No accessible voxel is available for particle initialization")
    selected = starts[rng.integers(0, len(starts), size=particle_count)]
    spacing = np.asarray((dz_um, dx_um, dx_um), dtype=float)
    positions = (selected.astype(float) + 0.5) * spacing
    escape_times = np.full(particle_count, np.nan)
    outcomes = np.zeros(particle_count, dtype=np.uint8)
    active = np.ones(particle_count, dtype=bool)
    step_sigma = np.sqrt(2.0 * diffusion_um2_s * dt_s)
    total_steps = int(np.ceil(max_time_s / dt_s))
    rejected_solid = rejected_lateral = 0
    for step in range(1, total_steps + 1):
        ids = np.flatnonzero(active)
        if not len(ids):
            break
        h = bounded_step_duration(step, dt_s, max_time_s)
        if h <= 0:
            break
        proposal = positions[ids] + rng.normal(
            0.0, np.sqrt(2.0*diffusion_um2_s*h), size=(len(ids), 3)
        )
        events, _ = trace_proposals(accessible, positions[ids], proposal, spacing)
        escaped = (events == ESCAPE_Z_MIN) | (events == ESCAPE_Z_MAX)
        if np.any(escaped):
            escaped_ids = ids[escaped]
            escape_times[escaped_ids] = min(step * dt_s, max_time_s)
            outcomes[escaped_ids] = events[escaped]
            active[escaped_ids] = False
        allowed = events == ACCEPT
        positions[ids[allowed]] = proposal[allowed]
        rejected_solid += int(np.count_nonzero(events == REJECT_SOLID))
        rejected_lateral += int(np.count_nonzero(events == REJECT_LATERAL))
    escaped_mask = np.isfinite(escape_times)
    observed = np.where(escaped_mask, escape_times, max_time_s)
    grid = np.linspace(0.0, max_time_s, 101)
    survival = empirical_survival(escape_times, grid, max_time_s)
    escape_fraction = float(np.mean(escaped_mask))
    retained_fraction = 1.0 - escape_fraction
    summary = {
        "particle_count": int(particle_count),
        "escaped_count": int(np.count_nonzero(escaped_mask)),
        "escaped_z_min_count": int(np.count_nonzero(outcomes == ESCAPE_Z_MIN)),
        "escaped_z_max_count": int(np.count_nonzero(outcomes == ESCAPE_Z_MAX)),
        "censored_count": int(np.count_nonzero(~escaped_mask)),
        "escape_fraction": escape_fraction,
        "retained_fraction_at_censor": retained_fraction,
        "restricted_mean_first_passage_time_s": float(np.mean(observed)),
        "MET_restricted_s": float(np.mean(observed)),
        "conditional_mean_escape_time_s": (
            float(np.mean(escape_times[escaped_mask])) if np.any(escaped_mask) else None
        ),
        "median_observed_or_censored_time_s": float(np.median(observed)),
        "right_censor_time_s": float(max_time_s),
        "initialization": initialization,
        "central_start_voxel_count": int(len(starts)),
        "solid_collision_rejections": rejected_solid,
        "lateral_boundary_rejections": rejected_lateral,
        "step_sigma_per_axis_um": float(step_sigma),
        "transport_kernel_version": TRANSPORT_KERNEL_VERSION,
        "collision_rule": COLLISION_RULE,
        "last_step_rule": "variance 2*D*min(dt, remaining observation time)",
        "exit_definition": "unobstructed first particle-center contact with either z face",
        "escape_time_resolution": "end of transport step; not exact Brownian hitting time",
        "survival_definition": "mean(censored | (observed escape time > t)), 0<=t<=Tc",
        "MET_definition": (
            "mean(min(T_escape, censor_time)); this is a restricted mean first-passage "
            "time, not an uncensored arithmetic mean"
        ),
        "tracer_injection_rate_cells_s": float(injection_rate_cells_s),
        "aggregation_rate_Ragg_cells_s": (
            float(injection_rate_cells_s * retained_fraction)
            if injection_rate_cells_s > 0.0
            else None
        ),
        "Ragg_definition": (
            "R_agg=J_in*(1-P_escape). It is reported only when a nonzero external "
            "injection flux J_in is supplied; the batch simulation alone does not define cells/s."
        ),
    }
    result = (summary, escape_times, np.column_stack((grid, survival)))
    return (*result, outcomes) if return_outcomes else result
