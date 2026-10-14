"""Rest-angle and archive-literal planar bending forces."""

from __future__ import annotations

import numpy as np

from .config import Program19Config
from .state import SimulationState, rms_vectors


def _perp(vector: np.ndarray) -> np.ndarray:
    return np.asarray((-vector[1], vector[0]), dtype=float)


def _wrapped_difference(value: float, reference: float) -> float:
    return float(np.arctan2(np.sin(value - reference), np.cos(value - reference)))


def bending_triplet(
    xi: np.ndarray,
    xj: np.ndarray,
    xk: np.ndarray,
    phi0: float,
    dual_L0_um: float,
    config: Program19Config,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u, v = xj - xi, xk - xj
    lu2, lv2 = float(np.dot(u, u)), float(np.dot(v, v))
    if lu2 <= 1e-18 or lv2 <= 1e-18:
        zero = np.zeros(2)
        return zero.copy(), zero.copy(), zero.copy()
    phi = float(np.arctan2(u[0] * v[1] - u[1] * v[0], np.dot(u, v)))
    delta_phi = _wrapped_difference(phi, phi0)
    if config.equation_profile == "codeocean_literal":
        lu, lv = np.sqrt(lu2), np.sqrt(lv2)
        magnitude = config.archive_bending_force_nN * delta_phi
        fi = magnitude * (-_perp(u) / lu)
        fk = -magnitude * (_perp(v) / lv)
        fj = -fi - fk
    else:
        k_theta_nN_um = config.bending_EI_nN_um2 / max(dual_L0_um, 1e-12)
        dphi_du = -_perp(u) / lu2
        dphi_dv = _perp(v) / lv2
        dE_du = k_theta_nN_um * delta_phi * dphi_du
        dE_dv = k_theta_nN_um * delta_phi * dphi_dv
        fi = dE_du
        fj = -dE_du + dE_dv
        fk = -dE_dv
    return fi, fj, fk


def bending_force_field(
    state: SimulationState, config: Program19Config
) -> tuple[np.ndarray, float]:
    force = np.zeros_like(state.nodes_um)
    for triplet, phi0, dual in zip(
        state.bending_triples, state.bending_phi0_rad, state.bending_dual_L0_um
    ):
        i, j, k = (int(value) for value in triplet)
        fi, fj, fk = bending_triplet(
            state.nodes_um[i], state.nodes_um[j], state.nodes_um[k], phi0, dual, config
        )
        force[i] += fi
        force[j] += fj
        force[k] += fk
    return force, rms_vectors(force)
