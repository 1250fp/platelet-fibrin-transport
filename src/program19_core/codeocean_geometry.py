"""Read an explicit Michael CodeOcean XML network and project it onto XY.

This importer is a comparison tool.  It preserves node/link identities and
fixed-node identifiers, but a 2-D projection cannot preserve the original 3-D
torsion angles or out-of-plane edge lengths.
"""

from __future__ import annotations

from pathlib import Path
import hashlib
import math
import xml.etree.ElementTree as ET

import numpy as np

from .config import Program19Config
from .state import SimulationState


def default_codeocean_xml(config: Program19Config) -> Path:
    name = "data_pdf.xml" if config.parameter_profile == "codeocean_data_pdf" else "data.xml"
    return Path(__file__).resolve().parents[1] / "reference_inputs" / name


def resolve_codeocean_xml(config: Program19Config) -> Path:
    path = Path(config.codeocean_xml_path).expanduser() if config.codeocean_xml_path else default_codeocean_xml(config)
    if not path.is_file():
        raise FileNotFoundError(
            f"CodeOcean XML not found: {path}. Set codeocean_xml_path or restore reference_inputs/."
        )
    return path.resolve()


def _numbers(element: ET.Element, count: int) -> list[float]:
    values = [float(value) for value in (element.text or "").split()]
    if len(values) < count:
        raise ValueError(f"Expected {count} values in <{element.tag}>, received {values}")
    return values[:count]


def _angle(xi: np.ndarray, xj: np.ndarray, xk: np.ndarray) -> float:
    u, v = xj - xi, xk - xj
    return float(np.arctan2(u[0] * v[1] - u[1] * v[0], np.dot(u, v)))


def _filopodia_count(rng: np.random.Generator, config: Program19Config) -> int:
    if config.filopodia_count_model == "fixed":
        return int(config.filopodia_fixed_count)
    if config.filopodia_count_model == "lognormal_clipped_round":
        return int(
            np.clip(
                round(
                    rng.lognormal(
                        config.filopodia_lognormal_mu,
                        config.filopodia_lognormal_sigma,
                    )
                ),
                config.filopodia_min,
                config.filopodia_max,
            )
        )
    for _ in range(100_000):
        value = float(rng.lognormal(config.filopodia_lognormal_mu, config.filopodia_lognormal_sigma))
        if config.filopodia_min <= value < config.filopodia_max + 1:
            return int(value)
    raise RuntimeError("Could not reproduce the bounded CodeOcean filopodia draw")


def _fixed_ids(root: ET.Element) -> list[int]:
    fixed = root.find("fixed")
    if fixed is None:
        return []
    values: list[int] = []
    for child in list(fixed):
        text = (child.text or "").strip()
        if text:
            values.append(int(text.split()[0]))
    if not list(fixed) and (fixed.text or "").strip():
        values.extend(int(value) for value in (fixed.text or "").split())
    return values


def build_codeocean_projected_state(config: Program19Config) -> SimulationState:
    path = resolve_codeocean_xml(config)
    root = ET.parse(path).getroot()
    node_elements = root.findall("./nodes/node")
    link_elements = root.findall("./links/link")
    if not node_elements or not link_elements:
        raise ValueError(f"{path} does not contain the expected nodes and links")
    nodes_xyz = np.asarray([_numbers(node, 3) for node in node_elements], dtype=float)
    nodes = nodes_xyz[:, :2].copy()
    center = np.mean(nodes, axis=0)
    nodes -= center
    platelets_xyz = np.asarray(
        [_numbers(platelet, 3) for platelet in root.findall("./plts/plt")], dtype=float
    ).reshape((-1, 3))
    platelets = platelets_xyz[:, :2] - center if len(platelets_xyz) else np.empty((0, 2))
    pairs = np.asarray(
        [[int(value) for value in (link.text or "").split()[:2]] for link in link_elements],
        dtype=int,
    ).reshape((-1, 2))
    if np.min(pairs) < 0 or np.max(pairs) >= len(nodes):
        raise ValueError("A CodeOcean link refers to a node outside the XML node list")
    delta = nodes[pairs[:, 1]] - nodes[pairs[:, 0]]
    lengths = np.linalg.norm(delta, axis=1)
    if np.any(lengths <= 1e-12):
        raise ValueError(
            "At least one 3-D CodeOcean link collapses in the XY projection; choose another projection."
        )
    incident: list[list[int]] = [[] for _ in nodes]
    adjacency: list[list[int]] = [[] for _ in nodes]
    node_fibers: list[set[int]] = [set() for _ in nodes]
    for edge_id, (a_raw, b_raw) in enumerate(pairs):
        a, b = int(a_raw), int(b_raw)
        incident[a].append(edge_id)
        incident[b].append(edge_id)
        adjacency[a].append(b)
        adjacency[b].append(a)
        # CodeOcean does not use this identity to restrict linking.  Unique
        # edge IDs are retained only for provenance and visualization.
        node_fibers[a].add(edge_id)
        node_fibers[b].add(edge_id)

    triples: list[tuple[int, int, int]] = []
    rest_angles: list[float] = []
    dual_lengths: list[float] = []
    for middle, neighbors in enumerate(adjacency):
        if len(neighbors) < 2:
            continue
        ordered = sorted(
            neighbors,
            key=lambda node: math.atan2(
                nodes[node, 1] - nodes[middle, 1], nodes[node, 0] - nodes[middle, 0]
            ),
        )
        neighbor_pairs = (
            [(ordered[0], ordered[1])]
            if len(ordered) == 2
            else list(zip(ordered, ordered[1:] + ordered[:1]))
        )
        for left, right in neighbor_pairs:
            triples.append((left, middle, right))
            rest_angles.append(_angle(nodes[left], nodes[middle], nodes[right]))
            dual_lengths.append(
                0.5
                * (
                    np.linalg.norm(nodes[middle] - nodes[left])
                    + np.linalg.norm(nodes[right] - nodes[middle])
                )
            )

    fixed = np.zeros(len(nodes), dtype=bool)
    for node_id in _fixed_ids(root):
        if 0 <= node_id < len(fixed):
            fixed[node_id] = True
        else:
            raise ValueError(f"Fixed node id {node_id} is outside the XML node list")
    rng = np.random.default_rng(np.random.SeedSequence([config.seed, 4019]))
    filopodia = [[None] * _filopodia_count(rng, config) for _ in platelets]
    return SimulationState(
        nodes_um=nodes,
        platelets_um=platelets,
        backbone_pairs=pairs,
        backbone_L0_um=lengths,
        bending_triples=np.asarray(triples, dtype=int).reshape((-1, 3)),
        bending_phi0_rad=np.asarray(rest_angles, dtype=float),
        bending_dual_L0_um=np.asarray(dual_lengths, dtype=float),
        node_fiber_ids=[frozenset(values) for values in node_fibers],
        incident_backbone=incident,
        original_adjacency=adjacency,
        fixed_nodes=fixed,
        filopodia=filopodia,
    )


def codeocean_source_metadata(config: Program19Config) -> dict:
    if config.geometry_profile != "codeocean_xml_xy":
        return {"used": False}
    path = resolve_codeocean_xml(config)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    root = ET.parse(path).getroot()
    settings = {
        child.tag: (child.text or "").strip()
        for child in root.findall("./settings/*")
    }
    return {
        "used": True,
        "path": str(path),
        "sha256": digest,
        "projection": "XY centered by the projected node centroid",
        "node_count": len(root.findall("./nodes/node")),
        "link_count": len(root.findall("./links/link")),
        "platelet_count": len(root.findall("./plts/plt")),
        "settings": settings,
        "warning": "A 2-D projection cannot preserve original 3-D torsion geometry.",
    }
