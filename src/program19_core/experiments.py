"""Exploratory, reproducible parameter sweeps for the revised minipaper.

These routines call the same solver as a normal run but write their results to
an explicitly labelled ``exploratory_experiments`` tree.  This keeps one-use
analyses available without presenting them as calibrated thesis evidence.
"""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

from .config import Program19Config
from .profiles import profile_audit
from .mechanical_runner import run_pipeline


DEFAULT_F0_SWEEP_NN = (0.0, 2.2, 4.4, 8.8)


def paired_replicate_base_seed(config: Program19Config, replicate: int) -> int:
    """Return a non-overlapping layer-seed block for one paired replicate."""
    if replicate < 0:
        raise ValueError("replicate must be nonnegative")
    return int(
        config.seed
        + replicate * config.mechanical_layers * config.layer_seed_stride
    )


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _aggregate(rows: list[dict], f0_values: Iterable[float]) -> list[dict]:
    metrics = (
        "geometric_porosity",
        "accessible_porosity_for_particle_centers",
        "escape_fraction",
        "retained_fraction_at_censor",
        "MET_restricted_s",
        "conditional_mean_escape_time_s",
        "aggregation_rate_Ragg_cells_s",
    )
    result: list[dict] = []
    for f0 in f0_values:
        selected = [row for row in rows if np.isclose(row["filopodia_f0_nN"], f0)]
        item: dict = {"filopodia_f0_nN": f0, "replicate_count": len(selected)}
        for metric in metrics:
            values = np.asarray(
                [row[metric] for row in selected if row.get(metric) is not None], dtype=float
            )
            item[f"{metric}_mean"] = float(np.mean(values)) if len(values) else None
            item[f"{metric}_std"] = (
                float(np.std(values, ddof=1)) if len(values) > 1 else 0.0 if len(values) else None
            )
        result.append(item)
    return result


def _plot(path: Path, aggregate: list[dict]) -> None:
    import matplotlib.pyplot as plt

    x = np.asarray([row["filopodia_f0_nN"] for row in aggregate], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1), constrained_layout=True)
    panels = (
        (
            "accessible_porosity_for_particle_centers",
            "Accessible porosity",
            "Fraction of analyzed region",
        ),
        ("MET_restricted_s", "Restricted MET", "Time (s)"),
        ("aggregation_rate_Ragg_cells_s", r"$R_{agg}$", "Cells/s"),
    )
    for axis, (metric, title, ylabel) in zip(axes, panels):
        mean = np.asarray(
            [np.nan if row[f"{metric}_mean"] is None else row[f"{metric}_mean"] for row in aggregate]
        )
        std = np.asarray(
            [0.0 if row[f"{metric}_std"] is None else row[f"{metric}_std"] for row in aggregate]
        )
        axis.errorbar(x, mean, yerr=std, marker="o", linewidth=1.8, capsize=3)
        axis.set_title(title)
        axis.set_xlabel(r"Filopodial baseline $f_0$ (nN)")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
        if metric == "aggregation_rate_Ragg_cells_s" and np.all(~np.isfinite(mean)):
            axis.text(
                0.5,
                0.5,
                "Undefined: specify an\ninjection flux $J_{in}$",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
    fig.suptitle("Program19.3 exploratory contractility sweep")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run_f0_sweep(
    base_config: Program19Config,
    f0_values_nN: Iterable[float] = DEFAULT_F0_SWEEP_NN,
    replicates: int = 1,
    progress_callback: Callable[[dict], None] | None = None,
) -> Path:
    """Run paired-seed f0 simulations and save raw plus aggregate metrics."""
    f0_values = tuple(float(value) for value in f0_values_nN)
    if not f0_values or any(value < 0.0 for value in f0_values):
        raise ValueError("f0_values_nN must contain nonnegative values")
    if replicates < 1:
        raise ValueError("replicates must be at least 1")
    root = (
        Path(base_config.output_root).expanduser()
        if base_config.output_root
        else Path(__file__).resolve().parents[1] / "results"
    )
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    experiment = root / "exploratory_experiments" / f"f0_sweep_{timestamp}"
    runs_root = experiment / "runs"
    runs_root.mkdir(parents=True, exist_ok=False)
    rows: list[dict] = []
    total = len(f0_values) * replicates
    completed = 0
    for replicate in range(replicates):
        paired_seed = paired_replicate_base_seed(base_config, replicate)
        for f0 in f0_values:
            values = base_config.to_dict()
            values.update(
                filopodia_f0_nN=f0,
                seed=paired_seed,
                run_transport=True,
                output_root=str(runs_root),
            )
            config = Program19Config.from_dict(values)

            def relay(payload: dict) -> None:
                if progress_callback:
                    progress_callback(
                        {
                            **payload,
                            "experiment": "f0_sweep",
                            "f0_nN": f0,
                            "replicate": replicate,
                            "run_index": completed,
                            "run_count": total,
                        }
                    )

            _, run_dir, status = run_pipeline(config, progress_callback=relay)
            metrics = json.loads((run_dir / "run_metrics.json").read_text(encoding="utf-8"))
            rows.append(
                {
                    "filopodia_f0_nN": f0,
                    "replicate": replicate,
                    "seed": paired_seed,
                    "status": status,
                    "run_directory": str(run_dir),
                    **{key: value for key, value in metrics.items() if key != "definitions"},
                }
            )
            completed += 1
    aggregate = _aggregate(rows, f0_values)
    _write_csv(experiment / "f0_sweep_raw.csv", rows)
    _write_csv(experiment / "f0_sweep_aggregate.csv", aggregate)
    _plot(experiment / "f0_sweep_porosity_MET_Ragg.png", aggregate)
    manifest = {
        "experiment": "exploratory f0 sweep",
        "intended_use": (
            "Minipaper revision and one-use sensitivity analysis. This label does not imply "
            "experimental calibration or thesis validation."
        ),
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "f0_values_nN": list(f0_values),
        "replicates": replicates,
        "paired_seed_rule": (
            "replicate_base_seed = base_seed + replicate*mechanical_layers*"
            "layer_seed_stride; each f0 reuses the same non-overlapping layer-seed block"
        ),
        "paired_replicate_base_seeds": [
            paired_replicate_base_seed(base_config, replicate)
            for replicate in range(replicates)
        ],
        "profiles": profile_audit(base_config),
        "metric_definitions": {
            "geometric_porosity": "void voxels / analyzed-region voxels",
            "accessible_porosity": (
                "particle-center-accessible voxels after EDT erosion / analyzed-region voxels"
            ),
            "MET_restricted_s": (
                "mean(min(T_escape,T_censor)); right-censored restricted mean first-passage time"
            ),
            "Ragg_cells_s": (
                "J_in*(1-P_escape); undefined if tracer_injection_rate_cells_s is zero"
            ),
        },
        "base_configuration": base_config.to_dict(),
    }
    (experiment / "experiment_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return experiment
