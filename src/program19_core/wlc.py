"""Audited Michael/Marko--Siggia stretch law and local tangent."""

from __future__ import annotations

import numpy as np

from .config import Program19Config
from .state import SimulationState, rms_vectors

K_B_J_PER_K = 1.380649e-23
J_PER_NN_UM = 1.0e-15


def kbt_nN_um(config: Program19Config) -> float:
    return K_B_J_PER_K * config.temperature_K / J_PER_NN_UM


def engineering_strain(distance_um, reference_um, config: Program19Config):
    distance = np.asarray(distance_um, dtype=float)
    reference = np.maximum(np.asarray(reference_um, dtype=float), 1e-12)
    return np.minimum((distance - reference) / reference, config.gamma_max)


def wlc_force_nN(distance_um, reference_um, config: Program19Config):
    """F(gamma), where gamma=(L-L0)/L0 and u=gamma/C.

    Compression is signed. Only the upper constitutive evaluation is capped;
    coordinates are never projected and displacements are never clipped.
    """
    gamma = engineering_strain(distance_um, reference_um, config)
    normalized = gamma / config.contour_multiplier_C
    denominator = np.maximum(1.0 - normalized, 1e-12)
    prefactor = (
        config.monomers_in_parallel_N
        * kbt_nN_um(config)
        / config.persistence_length_um
    )
    result = prefactor * (0.25 / denominator**2 - 0.25 + normalized)
    return float(result) if np.ndim(result) == 0 else result


def paper_stiffness_dF_dgamma_nN(strain, config: Program19Config):
    gamma = np.minimum(np.asarray(strain, dtype=float), config.gamma_max)
    C = config.contour_multiplier_C
    denominator = np.maximum(C - gamma, 1e-12)
    prefactor = (
        config.monomers_in_parallel_N
        * kbt_nN_um(config)
        / config.persistence_length_um
    )
    result = prefactor * (C**2 / (2.0 * denominator**3) + 1.0 / C)
    return float(result) if np.ndim(result) == 0 else result


def archive_literal_local_stiffness_nN(strain, config: Program19Config):
    """Literal audited CUDA stiffness convention, not an analytical derivative.

    The archive averages absolute strain and uses a first-power denominator in
    its second term. This profile exists for reproduction and sensitivity only.
    """
    gamma = np.minimum(np.abs(np.asarray(strain, dtype=float)), config.gamma_max)
    C = config.contour_multiplier_C
    denominator = np.maximum(C - gamma, 1e-12)
    prefactor = (
        config.monomers_in_parallel_N
        * kbt_nN_um(config)
        / config.persistence_length_um
    )
    result = prefactor * (1.0 / C + C**2 / (2.0 * denominator))
    return float(result) if np.ndim(result) == 0 else result


def pair_force_field(
    nodes: np.ndarray,
    pairs: np.ndarray,
    reference_lengths: np.ndarray,
    config: Program19Config,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    force = np.zeros_like(nodes)
    if len(pairs) == 0:
        return force, np.empty(0), np.empty(0)
    delta = nodes[pairs[:, 1]] - nodes[pairs[:, 0]]
    distance = np.linalg.norm(delta, axis=1)
    safe = np.maximum(distance, 1e-12)
    magnitude = wlc_force_nN(distance, reference_lengths, config)
    vector = magnitude[:, None] * delta / safe[:, None]
    np.add.at(force, pairs[:, 0], vector)
    np.add.at(force, pairs[:, 1], -vector)
    strain = engineering_strain(distance, reference_lengths, config)
    return force, strain, magnitude


def backbone_and_cohesion_forces(
    state: SimulationState, config: Program19Config
) -> tuple[np.ndarray, float, float]:
    backbone, strain, _ = pair_force_field(
        state.nodes_um, state.backbone_pairs, state.backbone_L0_um, config
    )
    total = backbone
    if state.cohesion_pairs:
        pairs = np.asarray(state.cohesion_pairs, dtype=int).reshape((-1, 2))
        cohesive, _, _ = pair_force_field(
            state.nodes_um, pairs, np.asarray(state.cohesion_L0_um), config
        )
        total = total + cohesive
    maximum = float(np.max(np.abs(strain))) if len(strain) else 0.0
    return total, rms_vectors(total), maximum


def local_node_stiffness_nN(
    state: SimulationState, node_index: int, config: Program19Config
) -> float:
    incident = state.incident_backbone[node_index]
    if not incident:
        return 0.0
    edge_ids = np.asarray(incident, dtype=int)
    pairs = state.backbone_pairs[edge_ids]
    distance = np.linalg.norm(
        state.nodes_um[pairs[:, 1]] - state.nodes_um[pairs[:, 0]], axis=1
    )
    strain = (distance - state.backbone_L0_um[edge_ids]) / np.maximum(
        state.backbone_L0_um[edge_ids], 1e-12
    )
    if config.equation_profile == "codeocean_literal":
        return float(archive_literal_local_stiffness_nN(np.mean(np.abs(strain)), config))
    return float(paper_stiffness_dF_dgamma_nN(np.mean(strain), config))
