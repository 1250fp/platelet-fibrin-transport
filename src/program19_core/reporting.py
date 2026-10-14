"""Deterministic report figures produced directly from Program19 outputs."""

from __future__ import annotations

from pathlib import Path
import matplotlib
import numpy as np

# Program19 writes figures to disk from GUI worker threads and spawned layer
# processes; it never embeds a Matplotlib window.  A non-interactive backend
# avoids Tk objects being created or destroyed outside Tk's main thread.
matplotlib.use("Agg")

from .config import Program19Config
from .state import SimulationState
from .transport import (
    PORE_BLIND_Z_CONNECTED,
    PORE_NO_Z_FACE_CONTACT,
    PORE_THROUGH_STACK,
)


def save_mechanics_figure(
    path: Path, state: SimulationState, config: Program19Config, status: str
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Circle

    figure, axis = plt.subplots(figsize=(7.2, 7.2), dpi=170)
    if len(state.backbone_pairs):
        axis.add_collection(
            LineCollection(
                state.nodes_um[state.backbone_pairs],
                colors="#2b4c6f",
                linewidths=0.8,
                alpha=0.95,
            )
        )
    if state.cohesion_pairs:
        pairs = np.asarray(state.cohesion_pairs, dtype=int)
        axis.add_collection(
            LineCollection(
                state.nodes_um[pairs], colors="#d95f02", linewidths=1.0, alpha=0.9
            )
        )
    for center in state.platelets_um:
        axis.add_patch(
            Circle(
                center,
                config.platelet_radius_um,
                facecolor="#ceb07a",
                edgecolor="#5b4630",
                linewidth=0.7,
            )
        )
    half = 0.53 * config.domain_size_um
    axis.set(xlim=(-half, half), ylim=(-half, half), xlabel="x (µm)", ylabel="y (µm)")
    axis.set_aspect("equal")
    axis.set_title(
        f"Program19.3 · {config.equation_profile} · {config.parameter_profile} · "
        f"t = {state.time_s:.4g} s · {status}"
    )
    axis.grid(alpha=0.12)
    figure.tight_layout()
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)


def save_region_figure(
    path: Path,
    solid: np.ndarray,
    region: np.ndarray,
    extent: np.ndarray,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    # 0 external, 1 void within analyzed region, 2 solid.
    classes = np.zeros_like(solid, dtype=np.uint8)
    classes[region] = 1
    classes[solid & region] = 2
    cmap = ListedColormap(("#e8edf0", "#49a9dc", "#d8c49b"))
    figure, axis = plt.subplots(figsize=(7.2, 7.2), dpi=180)
    axis.imshow(classes, origin="lower", extent=extent, cmap=cmap, vmin=0, vmax=2)
    axis.contour(
        region.astype(float),
        levels=(0.5,),
        origin="lower",
        extent=extent,
        colors=("#e39b00",),
        linewidths=2.5,
    )
    axis.set(xlabel="x (µm)", ylabel="y (µm)", aspect="equal")
    axis.set_title("Analyzed clot region and excluded exterior")
    figure.tight_layout()
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)


def save_stack_slices_figure(
    path: Path,
    solid_stack: np.ndarray,
    pore_classes: np.ndarray,
    extent: np.ndarray,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    indices = sorted({0, len(solid_stack) // 2, len(solid_stack) - 1})
    figure, axes = plt.subplots(2, len(indices), figsize=(4.5 * len(indices), 8.0), dpi=160)
    if len(indices) == 1:
        axes = np.asarray(axes).reshape((2, 1))
    pore_cmap = ListedColormap(("#081722", "#24d0d0", "#ffad32", "#b85cdd"))
    for column, layer in enumerate(indices):
        axes[0, column].imshow(
            solid_stack[layer], origin="lower", extent=extent, cmap="Greys", vmin=0, vmax=1
        )
        axes[0, column].set_title(f"Solid phase · layer {layer + 1}")
        axes[1, column].imshow(
            pore_classes[layer],
            origin="lower",
            extent=extent,
            cmap=pore_cmap,
            vmin=0,
            vmax=3,
        )
        axes[1, column].set_title(f"Classified void · layer {layer + 1}")
        for axis in axes[:, column]:
            axis.set_aspect("equal")
            axis.set(xlabel="x (µm)", ylabel="y (µm)")
    figure.legend(
        handles=(
            Patch(color="#24d0d0", label="Through-stack void"),
            Patch(color="#ffad32", label="Blind z-connected void"),
            Patch(color="#b85cdd", label="Internal component — no z-face contact"),
        ),
        loc="lower center",
        ncol=3,
        frameon=False,
    )
    figure.tight_layout(rect=(0, 0.07, 1, 1))
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)
