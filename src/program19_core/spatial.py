"""Small neighbor-search abstraction with a NumPy fallback."""

from __future__ import annotations

import math
import numpy as np

try:
    from scipy.spatial import cKDTree
except ImportError:  # pragma: no cover
    cKDTree = None


class NeighborIndex:
    def __init__(self, points: np.ndarray, cell_size: float):
        self.points = np.asarray(points, dtype=float)
        self.cell_size = max(float(cell_size), 1e-12)
        self.tree = cKDTree(self.points) if cKDTree is not None and len(self.points) else None
        self.cells: dict[tuple[int, int], list[int]] = {}
        if self.tree is None:
            for index, coordinate in enumerate(np.floor(self.points / self.cell_size).astype(int)):
                self.cells.setdefault((int(coordinate[0]), int(coordinate[1])), []).append(index)

    def query_ball_point(self, point: np.ndarray, radius: float) -> list[int]:
        if self.tree is not None:
            return [int(value) for value in self.tree.query_ball_point(point, float(radius))]
        center = np.floor(np.asarray(point) / self.cell_size).astype(int)
        span = int(math.ceil(float(radius) / self.cell_size))
        candidates: list[int] = []
        for dx in range(-span, span + 1):
            for dy in range(-span, span + 1):
                candidates.extend(self.cells.get((int(center[0] + dx), int(center[1] + dy)), []))
        if not candidates:
            return []
        ids = np.asarray(candidates, dtype=int)
        keep = np.sum((self.points[ids] - point) ** 2, axis=1) <= float(radius) ** 2
        return [int(value) for value in ids[keep]]

    def query_pairs(self, radius: float) -> set[tuple[int, int]]:
        if self.tree is not None:
            return {(int(a), int(b)) for a, b in self.tree.query_pairs(float(radius))}
        result: set[tuple[int, int]] = set()
        for index, point in enumerate(self.points):
            for neighbor in self.query_ball_point(point, radius):
                if neighbor > index:
                    result.add((index, neighbor))
        return result
