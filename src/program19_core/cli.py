"""Command-line interface and GUI dispatch for Program19.3."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Optional

from .config import (
    EQUATION_PROFILES,
    GEOMETRY_PROFILES,
    NUMERICAL_PROFILES,
    PARAMETER_PROFILES,
    Program19Config,
)
from .profiles import make_profiled_config, profiles_markdown
from .mechanical_runner import run_pipeline


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description="Program19.3 evolving mechanical-clot pipeline")
    result.add_argument("--headless", action="store_true", help="run without the GUI")
    result.add_argument("--config", type=str, help="JSON configuration file")
    result.add_argument("--duration", type=float)
    result.add_argument("--dt", type=float)
    result.add_argument("--seed", type=int)
    result.add_argument("--platelets", type=int)
    result.add_argument("--profile", dest="equation_profile", choices=EQUATION_PROFILES,
                        help="backward-compatible alias for --equation-profile")
    result.add_argument("--equation-profile", dest="equation_profile", choices=EQUATION_PROFILES)
    result.add_argument("--parameter-profile", choices=PARAMETER_PROFILES)
    result.add_argument("--numerical-profile", choices=NUMERICAL_PROFILES)
    result.add_argument("--geometry-profile", choices=GEOMETRY_PROFILES)
    result.add_argument("--codeocean-xml", type=str)
    result.add_argument("--output-root", type=str)
    result.add_argument("--layers", type=int, help="number of independent mechanical z layers")
    result.add_argument(
        "--workers", type=int,
        help="parallel layer processes (1 keeps serial execution)",
    )
    result.add_argument("--layer-seed-stride", type=int)
    result.add_argument("--no-intermediate-masks", action="store_true")
    result.add_argument("--transport", action="store_true")
    result.add_argument("--self-test", action="store_true")
    result.add_argument("--write-default-config", type=str)
    result.add_argument("--write-profile-guide", type=str)
    result.add_argument("--f0-sweep", type=str, help="comma-separated f0 values in nN")
    result.add_argument("--replicates", type=int, default=1)
    result.add_argument("--injection-rate", type=float)
    result.add_argument("--catalog-root", action="append", help="list reusable runs below this root")
    result.add_argument("--continue-run", type=str, help="saved mechanical run to continue in time")
    result.add_argument("--additional-time", type=float, help="additional mechanical time in seconds")
    result.add_argument("--extend-stack", type=str, help="saved stack to extend with new layers")
    result.add_argument("--target-layers", type=int)
    result.add_argument("--first-new-seed", type=int)
    result.add_argument("--evolving-transport-run", type=str, help="saved time-resolved run to analyze")
    result.add_argument(
        "--frozen-transport-run", type=str,
        help="saved time-resolved run to analyze at one frozen mechanical frame",
    )
    result.add_argument("--release-time", type=float, default=0.0)
    result.add_argument("--observation-time", type=float)
    result.add_argument("--particle-radius", type=float)
    result.add_argument("--diffusion", type=float)
    result.add_argument("--particles", type=int)
    result.add_argument("--transport-seed", type=int)
    result.add_argument("--recorded-particles", type=int, default=25)
    result.add_argument("--no-frozen-comparison", action="store_true")
    return result


def _self_test() -> None:
    import numpy as np
    from .bending import bending_triplet
    from .transport import classify_void_components
    from .wlc import paper_stiffness_dF_dgamma_nN, wlc_force_nN

    config = Program19Config(
        duration_s=0.002,
        dt_s=0.001,
        platelet_count=0,
        raster_resolution=32,
        mechanical_layers=3,
    )
    L0 = 0.3
    assert abs(wlc_force_nN(L0, L0, config)) < 1e-12
    assert wlc_force_nN(1.1 * L0, L0, config) > 0
    assert wlc_force_nN(0.9 * L0, L0, config) < 0
    gamma = 0.2
    eps = 1e-6
    numeric = (
        wlc_force_nN((1 + gamma + eps) * L0, L0, config)
        - wlc_force_nN((1 + gamma - eps) * L0, L0, config)
    ) / (2 * eps)
    assert np.isclose(numeric, paper_stiffness_dF_dgamma_nN(gamma, config), rtol=2e-5)
    xi, xj, xk = np.asarray((0.0, 0.0)), np.asarray((1.0, 0.0)), np.asarray((2.0, 0.0))
    forces = bending_triplet(xi, xj, xk, 0.0, 1.0, config)
    assert np.allclose(np.sum(forces, axis=0), 0.0)
    assert np.allclose(forces, 0.0)
    stack = np.ones((3, 16, 16), dtype=bool)
    stack[:, 5:9, 5:9] = False
    region = np.ones_like(stack)
    classes, summary = classify_void_components(stack, region)
    assert classes.shape == stack.shape and summary["void_voxels"] > 0
    print("Program19.3 self-test passed")


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = parser().parse_args(argv)
    if args.self_test:
        _self_test()
        return 0
    if args.write_profile_guide:
        Path(args.write_profile_guide).write_text(profiles_markdown(), encoding="utf-8")
        print(f"Profile guide written to {args.write_profile_guide}")
        return 0
    if args.write_default_config:
        make_profiled_config(
            equation_profile=args.equation_profile or "paper_consistent",
            parameter_profile=args.parameter_profile or "thesis_audited",
            numerical_profile=args.numerical_profile or "corrected_adaptive",
            geometry_profile=args.geometry_profile or "generated_2d",
        ).save_json(args.write_default_config)
        print(f"Default configuration written to {args.write_default_config}")
        return 0
    if args.catalog_root:
        from .catalog import scan_runs

        for record in scan_runs(args.catalog_root):
            print(
                f"{record.created_local}\t{record.status}\t{record.mechanical_time_s:g}s\t"
                f"{record.completed_layers}/{record.requested_layers} layers\t"
                f"f0={record.f0_nN:g}\t{record.path}"
            )
        return 0
    if args.continue_run:
        if args.additional_time is None:
            raise SystemExit("--continue-run requires --additional-time")
        from .continuation import continue_mechanical_run

        output_root = args.output_root or str(Path(__file__).resolve().parents[1] / "results")
        directory = continue_mechanical_run(
            args.continue_run, additional_time_s=args.additional_time, output_root=output_root,
            mechanical_workers=args.workers,
        )
        print(f"Mechanical continuation completed\nResults: {directory}")
        return 0
    if args.extend_stack:
        if args.target_layers is None or args.first_new_seed is None:
            raise SystemExit("--extend-stack requires --target-layers and --first-new-seed")
        from .stack_extension import extend_stack

        output_root = args.output_root or str(Path(__file__).resolve().parents[1] / "results")
        directory = extend_stack(
            args.extend_stack,
            target_layer_count=args.target_layers,
            first_new_seed=args.first_new_seed,
            output_root=output_root,
            mechanical_workers=args.workers,
        )
        print(f"Stack extension completed\nResults: {directory}")
        return 0
    if args.frozen_transport_run:
        from .evolving_transport import run_frozen_transport_analysis

        source_config = Program19Config.load_json(Path(args.frozen_transport_run) / "config.json")
        if args.observation_time is None:
            raise SystemExit("--frozen-transport-run requires --observation-time")
        output_root = args.output_root or str(
            Path(__file__).resolve().parents[1] / "results" / "transport_analyses"
        )
        directory = run_frozen_transport_analysis(
            args.frozen_transport_run,
            output_root=output_root,
            release_time_s=args.release_time,
            observation_time_s=args.observation_time,
            particle_radius_um=(
                args.particle_radius
                if args.particle_radius is not None else source_config.particle_radius_um
            ),
            diffusion_um2_s=(
                args.diffusion
                if args.diffusion is not None else source_config.tracer_diffusion_um2_s
            ),
            particle_count=(
                args.particles
                if args.particles is not None else source_config.transport_particles
            ),
            dt_s=source_config.transport_dt_s,
            seed=(
                args.transport_seed
                if args.transport_seed is not None else source_config.seed + 9403
            ),
            recorded_particles=args.recorded_particles,
        )
        print(f"Frozen transport completed\nResults: {directory}")
        return 0
    if args.evolving_transport_run:
        from .evolving_transport import run_evolving_transport_comparison

        source_config = Program19Config.load_json(Path(args.evolving_transport_run) / "config.json")
        if args.observation_time is None:
            raise SystemExit("--evolving-transport-run requires --observation-time")
        output_root = args.output_root or str(
            Path(__file__).resolve().parents[1] / "results" / "transport_analyses"
        )
        directory = run_evolving_transport_comparison(
            args.evolving_transport_run,
            output_root=output_root,
            release_time_s=args.release_time,
            observation_time_s=args.observation_time,
            particle_radius_um=args.particle_radius or source_config.particle_radius_um,
            diffusion_um2_s=args.diffusion or source_config.tracer_diffusion_um2_s,
            particle_count=args.particles or source_config.transport_particles,
            dt_s=source_config.transport_dt_s,
            seed=args.transport_seed if args.transport_seed is not None else source_config.seed + 9303,
            recorded_particles=args.recorded_particles,
            include_frozen=not args.no_frozen_comparison,
        )
        print(f"Evolving transport completed\nResults: {directory}")
        return 0
    if not args.headless:
        from .gui import Program19App
        Program19App().run()
        return 0
    config = (
        Program19Config.load_json(args.config)
        if args.config
        else make_profiled_config(
            equation_profile=args.equation_profile or "paper_consistent",
            parameter_profile=args.parameter_profile or "thesis_audited",
            numerical_profile=args.numerical_profile or "corrected_adaptive",
            geometry_profile=args.geometry_profile or "generated_2d",
        )
    )
    values = config.to_dict()
    overrides = {
        "duration_s": args.duration,
        "dt_s": args.dt,
        "seed": args.seed,
        "platelet_count": args.platelets,
        "output_root": args.output_root,
        "codeocean_xml_path": args.codeocean_xml,
        "tracer_injection_rate_cells_s": args.injection_rate,
        "mechanical_layers": args.layers,
        "mechanical_workers": args.workers,
        "layer_seed_stride": args.layer_seed_stride,
    }
    values.update({key: value for key, value in overrides.items() if value is not None})
    if args.no_intermediate_masks:
        values["save_intermediate_masks"] = False
    if args.transport:
        values["run_transport"] = True
    config = Program19Config.from_dict(values)

    if args.f0_sweep:
        from .experiments import run_f0_sweep

        f0_values = [float(value.strip()) for value in args.f0_sweep.split(",") if value.strip()]
        directory = run_f0_sweep(config, f0_values, args.replicates)
        print(f"Exploratory sweep completed\nResults: {directory}")
        return 0

    def show(payload: dict) -> None:
        if payload.get("stage") == "mechanics":
            print(
                f"\r{100 * payload['fraction']:6.2f}%  "
                f"layer={payload.get('layer_index', 0)}/{payload.get('layer_count', 0)}  "
                f"t={payload.get('simulated_time_s', 0):.4g} s  "
                f"links={payload.get('cohesion_links', 0)}  "
                f"substeps={payload.get('internal_substeps', 1)}",
                end="",
                flush=True,
            )
        elif payload.get("stage") == "complete":
            print()

    _, directory, status = run_pipeline(config, progress_callback=show)
    print(f"Status: {status}\nResults: {directory}")
    return 0
