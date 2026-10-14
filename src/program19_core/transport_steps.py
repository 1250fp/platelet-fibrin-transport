"""Shared finite-step transport geometry and administrative-censor reporting.

DDA visits the complete voxel supercover of each straight Euler proposal in
physical z,y,x coordinates. It prevents tunnelling through inaccessible cells;
it does not make a polygonal path an exact reflected Brownian path. Refused
proposals consume time. Temporal/spatial convergence must still be measured.

Traversal reference: Amanatides & Woo (1987), doi:10.2312/egtp.19871000.
"""
from __future__ import annotations

import numpy as np

TRANSPORT_KERNEL_VERSION = "dda_supercover_v1_20260907"
ACCEPT = 0
ESCAPE_Z_MIN = 1
ESCAPE_Z_MAX = 2
REJECT_SOLID = 3
REJECT_LATERAL = 4
COLLISION_RULE = (
    "ordered anisotropic voxel supercover of the proposal segment; first "
    "inaccessible or lateral contact rejects the entire proposal and consumes "
    "time; unobstructed first z-face contact absorbs; blocked edge/corner "
    "ties take precedence over escape"
)


def trace_proposals(accessible, positions, proposals, spacing):
    """Return (event_code, first_event_segment_fraction) for N proposals.

    Codes 0=accept, 1/2=lower/upper-z escape, 3=solid rejection, 4=lateral
    rejection. Fraction is 1 for an accepted segment. Escape fractions describe
    geometry only: they are NOT exact Brownian first-passage times.

    Inaccessible cells include the particle-radius exclusion already applied
    to the center mask. No second radius dilation is performed here. Contact
    with a blocked cell at an edge/corner is rejected (conservative supercover).
    """
    mask = np.asarray(accessible, dtype=bool)
    p = np.asarray(positions, dtype=float)
    q = np.asarray(proposals, dtype=float)
    spacing = np.asarray(spacing, dtype=float)
    if mask.ndim != 3 or any(s < 1 for s in mask.shape):
        raise ValueError("accessible must be a nonempty 3-D mask")
    if p.ndim != 2 or p.shape[1] != 3 or q.shape != p.shape:
        raise ValueError("positions and proposals must have matching (N,3) shapes")
    if spacing.shape != (3,) or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0):
        raise ValueError("spacing must contain positive finite (dz,dy,dx)")
    if not np.all(np.isfinite(p)) or not np.all(np.isfinite(q)):
        raise ValueError("positions and proposals must be finite")
    n = len(p)
    status = np.zeros(n, dtype=np.uint8)
    fraction = np.ones(n, dtype=float)
    if not n:
        return status, fraction
    start = p / spacing
    end = q / spacing
    delta = end - start
    shape = np.asarray(mask.shape)
    if np.any((start < 0) | (start >= shape)):
        raise ValueError("all starting centers must lie inside the grid")
    cells = np.floor(start).astype(np.int64)
    directions = np.sign(delta).astype(np.int64)

    def cell_events(c):
        result = np.zeros(len(c), dtype=np.uint8)
        lateral = np.any((c[:, 1:] < 0) | (c[:, 1:] >= shape[1:]), axis=1)
        result[c[:, 0] < 0] = ESCAPE_Z_MIN
        result[c[:, 0] >= shape[0]] = ESCAPE_Z_MAX
        result[lateral] = REJECT_LATERAL
        valid = np.flatnonzero(result == ACCEPT)
        if len(valid):
            v = c[valid]
            result[valid[~mask[v[:, 0], v[:, 1], v[:, 2]]]] = REJECT_SOLID
        return result

    # A segment exactly on a voxel plane touches cells on both sides. This is
    # rare with Gaussian increments, but essential for deterministic edge tests.
    static = (delta == 0) & (start == cells)
    static_codes = static @ np.asarray([1, 2, 4], dtype=np.int64)

    def covered_events(c, ids):
        result = cell_events(c)
        if np.any(static_codes[ids]):
            for code in range(1, 8):
                use = np.flatnonzero((static_codes[ids] & code) == code)
                if len(use):
                    offset = np.asarray([bool(code & (1 << axis)) for axis in range(3)])
                    result[use] = np.maximum(result[use], cell_events(c[use] - offset))
        return result

    ids = np.arange(n)
    status[:] = covered_events(cells, ids)
    fraction[status != ACCEPT] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        boundary = cells + (directions > 0)
        t_next = np.where(delta != 0, (boundary-start)/delta, np.inf)
        t_delta = np.where(delta != 0, 1/np.abs(delta), np.inf)
    remaining = np.flatnonzero((status == ACCEPT) & (np.min(t_next, axis=1) <= 1.0))
    iterations = 0
    while len(remaining):
        iterations += 1
        if iterations > int(shape.sum()) + 3:
            raise RuntimeError("voxel traversal failed to terminate")
        event_t = np.min(t_next[remaining], axis=1)
        # Treat floating-point simultaneous plane intersections as ties.
        tolerance = 32*np.finfo(float).eps*np.maximum(1.0, np.abs(event_t))
        crossing = np.abs(t_next[remaining] - event_t[:, None]) <= tolerance[:, None]
        candidate = cells[remaining] + directions[remaining]*crossing
        events = covered_events(candidate, remaining)
        codes = crossing @ np.asarray([1,2,4], dtype=np.int64)
        tied = np.flatnonzero((codes & (codes-1)) != 0)
        if len(tied):
            # The full diagonal neighbor alone misses cells touched at a corner.
            for code in range(1,8):
                use = tied[((codes[tied] & code) == code) & (codes[tied] != code)]
                if len(use):
                    axes = np.asarray([bool(code & (1 << a)) for a in range(3)])
                    cc = cells[remaining[use]] + directions[remaining[use]]*axes
                    events[use] = np.maximum(events[use], covered_events(cc, remaining[use]))
        stopped = events != ACCEPT
        stop_ids = remaining[stopped]
        status[stop_ids] = events[stopped]
        fraction[stop_ids] = np.clip(event_t[stopped], 0.0, 1.0)
        continued = remaining[~stopped]
        cells[continued] = candidate[~stopped]
        t_next[continued] = np.where(
            crossing[~stopped], t_next[continued]+t_delta[continued], t_next[continued]
        )
        remaining = continued[np.min(t_next[continued], axis=1) <= 1.0 + 32*np.finfo(float).eps]
    return status, fraction


def empirical_survival(escape_times, time_grid, censor_time):
    """Survival with one common administrative horizon and NaN censored times."""
    times = np.asarray(escape_times, dtype=float)
    grid = np.asarray(time_grid, dtype=float)
    if times.ndim != 1 or not len(times) or grid.ndim != 1:
        raise ValueError("nonempty event times and a one-dimensional grid are required")
    if not np.isfinite(censor_time) or censor_time < 0:
        raise ValueError("censor_time must be finite and nonnegative")
    if np.any(~np.isfinite(grid)) or np.any((grid < 0) | (grid > censor_time)):
        raise ValueError("survival grid must stay inside the observation horizon")
    escaped = np.isfinite(times)
    if np.any(np.isinf(times)) or np.any((times[escaped] < 0) | (times[escaped] > censor_time)):
        raise ValueError("escape times must be observed within the horizon or NaN")
    return np.asarray([np.mean((~escaped) | (times > t)) for t in grid])


def bounded_step_duration(step, dt_s, horizon_s):
    """Duration of a one-based step, shortening the last Gaussian increment."""
    return min(dt_s, max(0.0, horizon_s - (step-1)*dt_s))
