"""One-way, inlet-to-outlet transport through a saved mechanical stack.

This module is intentionally separate from the original center-release first
passage experiment.  It starts tracers at the z-min face, counts transmission
only at the z-max face, and keeps upstream return distinct from finite-time
retention.  The distinction is required before an imposed injection flux can
be converted into rates.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np


INLET_BOUNDARIES = ("return", "reflecting")


def through_clot_first_passage(
    accessible: np.ndarray,
    dx_um: float,
    dz_um: float,
    diffusion_um2_s: float,
    particle_count: int,
    dt_s: float,
    max_time_s: float,
    rng: np.random.Generator,
    injection_rates_cells_s: Iterable[float] = (1.0,),
    inlet_boundary: str = "return",
    inlet_disk_fraction: float = 0.20,
) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    """Simulate transport from the z-min inlet to the z-max outlet.

    ``return`` treats a crossing of z<0 as upstream return.  ``reflecting``
    conditions on entry by reflecting that proposal at z=0.  In both cases a
    crossing of z>=Lz is transmission, x/y exits and solid collisions are
    rejected, and active particles at ``max_time_s`` are right-censored.

    Returns ``(summary, terminal_times, outcome_codes, survival_curve)``.
    Outcome codes are 0=retained/censored, 1=transmitted, 2=returned.
    """
    accessible = np.asarray(accessible, dtype=bool)
    if accessible.ndim != 3 or any(size < 1 for size in accessible.shape):
        raise ValueError("accessible must be a non-empty (z,y,x) Boolean stack")
    if inlet_boundary not in INLET_BOUNDARIES:
        raise ValueError(f"inlet_boundary must be one of {INLET_BOUNDARIES}")
    for name, value in (
        ("dx_um", dx_um),
        ("dz_um", dz_um),
        ("diffusion_um2_s", diffusion_um2_s),
        ("dt_s", dt_s),
        ("max_time_s", max_time_s),
    ):
        if not np.isfinite(value) or float(value) <= 0.0:
            raise ValueError(f"{name} must be positive")
    if particle_count < 1:
        raise ValueError("particle_count must be at least 1")
    if not 0.0 < inlet_disk_fraction <= 0.5:
        raise ValueError("inlet_disk_fraction must be in (0, 0.5]")

    injection_rates = tuple(float(value) for value in injection_rates_cells_s)
    if not injection_rates or any(not np.isfinite(value) or value < 0.0 for value in injection_rates):
        raise ValueError("injection rates must contain finite nonnegative values")

    nz, ny, nx = accessible.shape
    yy, xx = np.indices((ny, nx))
    central_xy = (
        (xx - 0.5 * (nx - 1)) ** 2 + (yy - 0.5 * (ny - 1)) ** 2
        <= (inlet_disk_fraction * min(nx, ny)) ** 2
    )
    starts = np.argwhere(accessible[0] & central_xy)
    start_region = "central inlet disk"
    if not len(starts):
        starts = np.argwhere(accessible[0])
        start_region = "entire accessible inlet face (central disk was empty)"
    if not len(starts):
        raise RuntimeError("No accessible particle-center position exists on the inlet z face")

    selected = starts[rng.integers(0, len(starts), size=particle_count)]
    positions = np.empty((particle_count, 3), dtype=float)
    positions[:, 0] = 0.5 * dz_um
    positions[:, 1] = (selected[:, 0].astype(float) + 0.5) * dx_um
    positions[:, 2] = (selected[:, 1].astype(float) + 0.5) * dx_um

    spacing = np.asarray((dz_um, dx_um, dx_um), dtype=float)
    upper = np.asarray((nz * dz_um, ny * dx_um, nx * dx_um), dtype=float)
    step_sigma = float(np.sqrt(2.0 * diffusion_um2_s * dt_s))
    total_steps = int(np.ceil(max_time_s / dt_s))
    terminal_times = np.full(particle_count, np.nan)
    outcomes = np.zeros(particle_count, dtype=np.uint8)
    active = np.ones(particle_count, dtype=bool)

    for step in range(1, total_steps + 1):
        ids = np.flatnonzero(active)
        if not len(ids):
            break
        proposal = positions[ids] + rng.normal(0.0, step_sigma, size=(len(ids), 3))

        if inlet_boundary == "reflecting":
            below = proposal[:, 0] < 0.0
            proposal[below, 0] = np.abs(proposal[below, 0])
            returned = np.zeros(len(ids), dtype=bool)
        else:
            returned = proposal[:, 0] < 0.0

        transmitted = (~returned) & (proposal[:, 0] >= upper[0])
        terminal = returned | transmitted
        if np.any(terminal):
            terminal_ids = ids[terminal]
            terminal_times[terminal_ids] = min(step * dt_s, max_time_s)
            outcomes[ids[returned]] = 2
            outcomes[ids[transmitted]] = 1
            active[terminal_ids] = False

        staying_ids = ids[~terminal]
        staying = proposal[~terminal]
        if len(staying_ids):
            inside = (
                (staying[:, 0] >= 0.0)
                & (staying[:, 0] < upper[0])
                & (staying[:, 1] >= 0.0)
                & (staying[:, 1] < upper[1])
                & (staying[:, 2] >= 0.0)
                & (staying[:, 2] < upper[2])
            )
            candidates = np.floor(staying / spacing).astype(int)
            allowed = np.zeros(len(staying_ids), dtype=bool)
            valid = np.flatnonzero(inside)
            if len(valid):
                c = candidates[valid]
                allowed[valid] = accessible[c[:, 0], c[:, 1], c[:, 2]]
            positions[staying_ids[allowed]] = staying[allowed]

    transmitted_mask = outcomes == 1
    returned_mask = outcomes == 2
    retained_mask = outcomes == 0
    observed = np.where(np.isfinite(terminal_times), terminal_times, max_time_s)
    grid = np.linspace(0.0, max_time_s, 101)
    survival = np.asarray([np.mean(observed > time) for time in grid])

    p_transmitted = float(np.mean(transmitted_mask))
    p_returned = float(np.mean(returned_mask))
    p_retained = float(np.mean(retained_mask))
    restricted_residence = float(np.mean(observed))
    flux_scenarios = []
    for rate in injection_rates:
        flux_scenarios.append(
            {
                "J_in_cells_s": rate,
                "transmission_rate_cells_s": rate * p_transmitted,
                "upstream_return_rate_cells_s": rate * p_returned,
                "Ragg_retained_at_censor_cells_s": rate * p_retained,
                "legacy_nontransmission_rate_cells_s": rate * (1.0 - p_transmitted),
                "mean_resident_particles": rate * restricted_residence,
            }
        )

    summary = {
        "particle_count": int(particle_count),
        "transmitted_count": int(np.count_nonzero(transmitted_mask)),
        "returned_count": int(np.count_nonzero(returned_mask)),
        "retained_count_at_censor": int(np.count_nonzero(retained_mask)),
        "transmission_probability": p_transmitted,
        "upstream_return_probability": p_returned,
        "retained_probability_at_censor": p_retained,
        "probability_mass_check": p_transmitted + p_returned + p_retained,
        "restricted_mean_residence_time_s": restricted_residence,
        "conditional_mean_transmission_time_s": (
            float(np.mean(terminal_times[transmitted_mask])) if np.any(transmitted_mask) else None
        ),
        "conditional_mean_return_time_s": (
            float(np.mean(terminal_times[returned_mask])) if np.any(returned_mask) else None
        ),
        "median_terminal_or_censored_time_s": float(np.median(observed)),
        "right_censor_time_s": float(max_time_s),
        "diffusion_um2_s": float(diffusion_um2_s),
        "dt_s": float(dt_s),
        "step_sigma_per_axis_um": step_sigma,
        "stack_thickness_um": float(nz * dz_um),
        "inlet_boundary": inlet_boundary,
        "initialization": start_region,
        "outlet_definition": "crossing z-max face",
        "side_boundary_definition": "x/y proposals outside the domain are rejected",
        "outcome_codes": {"0": "retained_at_censor", "1": "transmitted", "2": "returned"},
        "flux_scenarios": flux_scenarios,
        "rate_interpretation": {
            "Ragg_retained_at_censor_cells_s": (
                "J_in*P(still inside at censor time); finite-horizon retention, not proven binding"
            ),
            "legacy_nontransmission_rate_cells_s": (
                "J_in*(1-P_transmitted); includes upstream returns when inlet_boundary=return"
            ),
            "mean_resident_particles": "J_in*restricted mean residence time",
        },
    }
    return summary, terminal_times, outcomes, np.column_stack((grid, survival))

