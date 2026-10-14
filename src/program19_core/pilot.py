"""Resume an existing f0 sweep to a paired five-replicate transport pilot."""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

from .config import Program19Config
from .experiments import paired_replicate_base_seed
from .mechanical_runner import run_pipeline
from .through_clot import through_clot_first_passage


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


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
        "transmission_probability",
        "upstream_return_probability",
        "retained_probability_at_censor",
        "restricted_mean_residence_time_s",
        "conditional_mean_transmission_time_s",
    )
    result: list[dict] = []
    for f0 in f0_values:
        selected = [row for row in rows if np.isclose(float(row["filopodia_f0_nN"]), f0)]
        item: dict = {"filopodia_f0_nN": f0, "replicate_count": len(selected)}
        for metric in metrics:
            values = np.asarray(
                [float(row[metric]) for row in selected if row.get(metric) not in (None, "")],
                dtype=float,
            )
            item[f"{metric}_mean"] = float(np.mean(values)) if len(values) else None
            item[f"{metric}_std"] = (
                float(np.std(values, ddof=1)) if len(values) > 1 else 0.0 if len(values) else None
            )
        result.append(item)
    return result


def _aggregate_flux(flux_rows: list[dict]) -> list[dict]:
    metrics = (
        "transmission_rate_cells_s",
        "upstream_return_rate_cells_s",
        "Ragg_retained_at_censor_cells_s",
        "legacy_nontransmission_rate_cells_s",
        "mean_resident_particles",
    )
    groups = sorted(
        {
            (float(row["filopodia_f0_nN"]), float(row["J_in_cells_s"]))
            for row in flux_rows
        }
    )
    result: list[dict] = []
    for f0, rate in groups:
        selected = [
            row
            for row in flux_rows
            if np.isclose(float(row["filopodia_f0_nN"]), f0)
            and np.isclose(float(row["J_in_cells_s"]), rate)
        ]
        item: dict = {
            "filopodia_f0_nN": f0,
            "J_in_cells_s": rate,
            "replicate_count": len(selected),
        }
        for metric in metrics:
            values = np.asarray([float(row[metric]) for row in selected], dtype=float)
            item[f"{metric}_mean"] = float(np.mean(values))
            item[f"{metric}_std"] = (
                float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
            )
        result.append(item)
    return result


def _load_stack_transport_inputs(run_directory: Path) -> tuple[np.ndarray, float, float]:
    with np.load(run_directory / "pore_analysis.npz") as pore:
        accessible = np.asarray(pore["accessible_centers"], dtype=bool)
    with np.load(run_directory / "mechanical_spatial_stack.npz") as stack:
        if str(stack["axis_order"].item()) != "z,y,x":
            raise RuntimeError(f"Unexpected axis order in {run_directory}")
        dx_um = float(stack["dx_um"].item())
        dz_um = float(stack["dz_um"].item())
    return accessible, dx_um, dz_um


def _f0_tag(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def _plot_pilot(path: Path, aggregate: list[dict]) -> None:
    import matplotlib.pyplot as plt

    x = np.asarray([row["filopodia_f0_nN"] for row in aggregate], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1), constrained_layout=True)
    panels = (
        ("transmission_probability", "Outlet transmission", "Probability"),
        ("retained_probability_at_censor", "Retained at censor", "Probability"),
        ("restricted_mean_residence_time_s", "Restricted residence time", "Time (s)"),
    )
    for axis, (metric, title, ylabel) in zip(axes, panels):
        mean = np.asarray([row[f"{metric}_mean"] for row in aggregate], dtype=float)
        std = np.asarray([row[f"{metric}_std"] for row in aggregate], dtype=float)
        axis.errorbar(x, mean, yerr=std, marker="o", linewidth=1.8, capsize=3)
        axis.set_title(title)
        axis.set_xlabel(r"Filopodial baseline $f_0$ (nN)")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
    fig.suptitle("Program19.2 inlet-to-outlet paired pilot")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run_through_clot_pilot(
    existing_sweep: str | Path,
    target_replicates: int = 5,
    injection_rates_cells_s: Iterable[float] = (1.0, 10.0, 100.0),
    inlet_boundary: str = "return",
    particle_count: int | None = None,
    resume_pilot: str | Path | None = None,
    progress_callback: Callable[[dict], None] | None = None,
) -> Path:
    """Reuse old masks, add only missing mechanics, and run one-way transport."""
    existing_sweep = Path(existing_sweep).expanduser().resolve()
    manifest_path = existing_sweep / "experiment_manifest.json"
    raw_path = existing_sweep / "f0_sweep_raw.csv"
    if not manifest_path.exists() or not raw_path.exists():
        raise FileNotFoundError("Select a Program19.2 f0_sweep folder containing its manifest and raw CSV")
    if target_replicates < 1:
        raise ValueError("target_replicates must be at least 1")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    base_config = Program19Config.from_dict(manifest["base_configuration"])
    f0_values = tuple(float(value) for value in manifest["f0_values_nN"])
    injection_rates = tuple(float(value) for value in injection_rates_cells_s)
    particles = int(particle_count or base_config.transport_particles)
    old_rows = _read_csv(raw_path)
    existing_runs = {
        (int(row["replicate"]), float(row["filopodia_f0_nN"])): Path(row["run_directory"])
        for row in old_rows
        if row.get("status") == "completed"
    }

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pilot_manifest = {
        "experiment": "Program19.2 paired inlet-to-outlet transport pilot",
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "source_f0_sweep": str(existing_sweep),
        "target_paired_replicates": target_replicates,
        "f0_values_nN": list(f0_values),
        "paired_seed_rule": (
            "base_seed + replicate*mechanical_layers*layer_seed_stride; same seed block at every f0"
        ),
        "paired_replicate_base_seeds": [
            paired_replicate_base_seed(base_config, replicate)
            for replicate in range(target_replicates)
        ],
        "existing_mechanical_replicates_reused": sorted(
            {replicate for replicate, _ in existing_runs if replicate < target_replicates}
        ),
        "new_mechanical_replicates_requested": [
            replicate
            for replicate in range(target_replicates)
            if not all((replicate, f0) in existing_runs for f0 in f0_values)
        ],
        "transport": {
            "initialization": "accessible particle centers on z-min inlet central disk",
            "outlet": "z-max only",
            "inlet_boundary": inlet_boundary,
            "particle_count_per_stack": particles,
            "particle_radius_um": base_config.particle_radius_um,
            "diffusion_um2_s": base_config.tracer_diffusion_um2_s,
            "dt_s": base_config.transport_dt_s,
            "censor_time_s": base_config.transport_max_time_s,
            "J_in_scenarios_cells_s": list(injection_rates),
        },
        "rate_warning": (
            "Ragg_retained_at_censor=J_in*P(still inside at censor). The legacy "
            "J_in*(1-P_transmitted) also includes upstream returns when inlet_boundary=return."
        ),
        "base_configuration": base_config.to_dict(),
    }
    if resume_pilot:
        pilot_root = Path(resume_pilot).expanduser().resolve()
        saved_manifest_path = pilot_root / "pilot_manifest.json"
        if not saved_manifest_path.exists():
            raise FileNotFoundError("The resume folder does not contain pilot_manifest.json")
        saved = json.loads(saved_manifest_path.read_text(encoding="utf-8"))
        if Path(saved["source_f0_sweep"]).resolve() != existing_sweep:
            raise ValueError("The selected resume pilot belongs to a different f0 sweep")
        if int(saved["target_paired_replicates"]) != target_replicates:
            raise ValueError("Target replicate count differs from the saved pilot")
        saved_transport = saved["transport"]
        if saved_transport["inlet_boundary"] != inlet_boundary:
            raise ValueError("Inlet boundary differs from the saved pilot")
        if int(saved_transport["particle_count_per_stack"]) != particles:
            raise ValueError("Particle count differs from the saved pilot")
        if tuple(float(value) for value in saved_transport["J_in_scenarios_cells_s"]) != injection_rates:
            raise ValueError("J_in scenarios differ from the saved pilot")
    else:
        pilot_root = existing_sweep.parent / f"through_clot_pilot_{target_replicates}rep_{timestamp}"
        pilot_root.mkdir(parents=True, exist_ok=False)
        (pilot_root / "pilot_manifest.json").write_text(
            json.dumps(pilot_manifest, indent=2, ensure_ascii=False), encoding="utf-8"
        )

    new_runs_root = pilot_root / "new_mechanical_runs"
    transport_root = pilot_root / "through_clot_transport"
    new_runs_root.mkdir(parents=True, exist_ok=True)
    transport_root.mkdir(parents=True, exist_ok=True)
    raw_output = pilot_root / "through_clot_raw.csv"
    flux_output = pilot_root / "through_clot_flux_scenarios_raw.csv"
    rows: list[dict] = _read_csv(raw_output) if raw_output.exists() else []
    flux_rows: list[dict] = _read_csv(flux_output) if flux_output.exists() else []
    completed_keys = {
        (int(row["replicate"]), float(row["filopodia_f0_nN"])) for row in rows
    }

    # A run may have completed immediately before an interruption but before its
    # transport row was written. Recover it by its saved config instead of
    # recomputing five expensive mechanical layers.
    for candidate in new_runs_root.iterdir():
        config_path = candidate / "config.json"
        status_path = candidate / "status.json"
        if not config_path.exists() or not status_path.exists():
            continue
        status_data = json.loads(status_path.read_text(encoding="utf-8"))
        if status_data.get("status") != "completed":
            continue
        config_data = json.loads(config_path.read_text(encoding="utf-8"))
        candidate_seed = int(config_data["seed"])
        block = base_config.mechanical_layers * base_config.layer_seed_stride
        offset = candidate_seed - base_config.seed
        if offset < 0 or offset % block:
            continue
        replicate = offset // block
        f0 = float(config_data["filopodia_f0_nN"])
        existing_runs.setdefault((replicate, f0), candidate)
    total = target_replicates * len(f0_values)
    completed = 0
    for replicate in range(target_replicates):
        seed = paired_replicate_base_seed(base_config, replicate)
        for f0 in f0_values:
            key = (replicate, f0)
            if key in completed_keys:
                completed += 1
                if progress_callback:
                    progress_callback(
                        {
                            "stage": "pilot_stack_complete",
                            "experiment": "through_clot_pilot",
                            "replicate": replicate,
                            "f0_nN": f0,
                            "stack_index": completed,
                            "stack_count": total,
                            "reused_completed_transport": True,
                        }
                    )
                continue
            if key in existing_runs:
                run_directory = existing_runs[key]
                mechanics_source = "reused_existing_f0_sweep"
            else:
                values = base_config.to_dict()
                values.update(
                    seed=seed,
                    filopodia_f0_nN=f0,
                    run_transport=False,
                    output_root=str(new_runs_root),
                )
                config = Program19Config.from_dict(values)

                def relay(payload: dict) -> None:
                    if progress_callback:
                        progress_callback(
                            {
                                **payload,
                                "experiment": "through_clot_pilot",
                                "replicate": replicate,
                                "f0_nN": f0,
                                "stack_index": completed,
                                "stack_count": total,
                            }
                        )

                _, run_directory, status = run_pipeline(config, progress_callback=relay)
                if status != "completed":
                    raise RuntimeError(f"Mechanical run did not complete: {run_directory}")
                mechanics_source = "new_missing_replicate"

            accessible, dx_um, dz_um = _load_stack_transport_inputs(run_directory)
            rng = np.random.default_rng(
                np.random.SeedSequence([seed, 4001, int(round(1000.0 * f0))])
            )
            summary, terminal_times, outcomes, survival = through_clot_first_passage(
                accessible=accessible,
                dx_um=dx_um,
                dz_um=dz_um,
                diffusion_um2_s=base_config.tracer_diffusion_um2_s,
                particle_count=particles,
                dt_s=base_config.transport_dt_s,
                max_time_s=base_config.transport_max_time_s,
                rng=rng,
                injection_rates_cells_s=injection_rates,
                inlet_boundary=inlet_boundary,
            )
            destination = transport_root / f"replicate_{replicate:03d}_f0_{_f0_tag(f0)}"
            destination.mkdir(parents=True, exist_ok=False)
            np.savez_compressed(
                destination / "through_clot_first_passage.npz",
                terminal_times_s=terminal_times,
                outcome_codes=outcomes,
                survival_time_s=survival[:, 0],
                survival_probability=survival[:, 1],
            )
            (destination / "through_clot_summary.json").write_text(
                json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
            )
            pore_summary = json.loads((run_directory / "pore_summary.json").read_text(encoding="utf-8"))
            row = {
                "filopodia_f0_nN": f0,
                "replicate": replicate,
                "seed": seed,
                "mechanics_source": mechanics_source,
                "run_directory": str(run_directory),
                "transport_directory": str(destination),
                "geometric_porosity": pore_summary.get("geometric_porosity"),
                "accessible_porosity_for_particle_centers": pore_summary.get(
                    "accessible_porosity_for_particle_centers"
                ),
                **{key: value for key, value in summary.items() if not isinstance(value, (dict, list))},
            }
            rows.append(row)
            for scenario in summary["flux_scenarios"]:
                flux_rows.append(
                    {
                        "filopodia_f0_nN": f0,
                        "replicate": replicate,
                        "seed": seed,
                        "inlet_boundary": inlet_boundary,
                        **scenario,
                    }
                )
            _write_csv(raw_output, rows)
            _write_csv(flux_output, flux_rows)
            completed += 1
            if progress_callback:
                progress_callback(
                    {
                        "stage": "pilot_stack_complete",
                        "experiment": "through_clot_pilot",
                        "replicate": replicate,
                        "f0_nN": f0,
                        "stack_index": completed,
                        "stack_count": total,
                        "run_directory": str(run_directory),
                    }
                )

    aggregate = _aggregate(rows, f0_values)
    _write_csv(pilot_root / "through_clot_aggregate.csv", aggregate)
    _write_csv(
        pilot_root / "through_clot_flux_scenarios_aggregate.csv",
        _aggregate_flux(flux_rows),
    )
    _plot_pilot(pilot_root / "through_clot_pilot_summary.png", aggregate)
    (pilot_root / "status.json").write_text(
        json.dumps({"status": "completed", "completed_stacks": completed}, indent=2),
        encoding="utf-8",
    )
    return pilot_root
