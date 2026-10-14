"""Validated configuration for Program19.3.

All mechanics use micrometres (um), nanonewtons (nN), and seconds (s).
Four independent profile axes prevent a source-reproduction choice from being
silently mixed with thesis parameters or corrected numerics:

``equation_profile`` selects constitutive equations, ``parameter_profile``
records the numerical parameter set, ``numerical_profile`` selects the time
integration/noise convention, and ``geometry_profile`` selects generated 2-D
geometry or an explicit XY projection of a Michael CodeOcean XML network.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any
import json
import math


EQUATION_PROFILES = ("paper_consistent", "codeocean_literal")
PARAMETER_PROFILES = (
    "thesis_audited",
    "minipaper_revision",
    "codeocean_data_xml",
    "codeocean_data_pdf",
    "custom",
)
NUMERICAL_PROFILES = ("corrected_adaptive", "codeocean_fixed_step")
GEOMETRY_PROFILES = ("generated_2d", "codeocean_xml_xy")
# Backward-compatible import name for scripts written for Program19.
PHYSICS_PROFILES = EQUATION_PROFILES
MOBILITY_MODELS = ("effective_network", "free_stokes")
THERMAL_NOISE_MODELS = ("gaussian_independent", "codeocean_uniform_shared")
FILOPODIA_COUNT_MODELS = (
    "lognormal_clipped_round",
    "codeocean_lognormal_rejection",
    "fixed",
)
FILOPODIA_SELECTION_RULES = ("nearest_free", "codeocean_random")
FILOPODIA_FORCE_DISTRIBUTIONS = ("full_per_arm", "divide_by_attached")
FILOPODIA_HANDOVER_RULES = ("global_reacquire", "graph_neighbor")
TRAJECTORY_DETAILS = ("positions", "topology")


@dataclass
class Program19Config:
    # Run controls
    duration_s: float = 1.0
    dt_s: float = 0.001
    output_interval_s: float = 0.10
    preview_interval_s: float = 0.10
    seed: int = 2035
    equation_profile: str = "paper_consistent"
    parameter_profile: str = "thesis_audited"
    numerical_profile: str = "corrected_adaptive"
    geometry_profile: str = "generated_2d"
    codeocean_xml_path: str = ""
    output_root: str = ""
    save_trajectory: bool = True
    trajectory_detail: str = "topology"

    # Initial planar clot geometry
    domain_size_um: float = 47.6
    effective_thickness_um: float = 1.0
    rho_degree3_um3: float = 0.02000
    rho_degree4_um3: float = 0.00379
    main_node_min_spacing_um: float = 1.0
    segment_target_um: float = 0.30
    platelet_count: int = 45
    platelet_near_node_fraction: float = 0.50

    # WLC: gamma=(L-L0)/L0 and u=gamma/C
    temperature_K: float = 300.0
    persistence_length_um: float = 0.001074033502990
    contour_multiplier_C: float = 2.351243606987286
    monomers_in_parallel_N: float = 50.0
    gamma_max: float = 2.0

    # Bending and cohesion
    bending_EI_nN_um2: float = 0.002262593420685
    archive_bending_force_nN: float = 0.002262593420685
    cohesion_distance_um: float = 0.10
    cohesion_search_stride: int = 10
    cohesion_max_new_per_node_per_search: int = 4
    cohesion_max_total_neighbors: int = 4
    cohesion_exclude_same_fiber: bool = True
    cohesion_degree_includes_backbone: bool = False

    # Platelets and active Michael contact branches
    platelet_radius_um: float = 1.13
    platelet_core_fraction: float = 0.90
    fiber_radius_um: float = 0.05
    platelet_force_radius_um: float = 2.18
    platelet_fibrin_force_nN: float = 4.0
    platelet_platelet_force_parameter: float = 4.0
    filopodia_lognormal_mu: float = 1.85
    filopodia_lognormal_sigma: float = 0.27
    filopodia_min: int = 2
    filopodia_max: int = 13
    filopodia_count_model: str = "lognormal_clipped_round"
    filopodia_fixed_count: int = 3
    filopodia_selection_rule: str = "nearest_free"
    filopodia_force_distribution: str = "full_per_arm"
    filopodia_global_node_exclusivity: bool = True
    filopodia_handover_rule: str = "global_reacquire"
    filopodia_f0_nN: float = 2.233
    filopodia_a: float = 61.0
    filopodia_b: float = 3.582
    filopodia_c: float = 131.8

    # Overdamped motion. The effective values are interpreted as nN*s/um.
    mobility_model: str = "effective_network"
    viscosity_Pa_s: float = 0.004
    effective_node_drag_nN_s_per_um: float = 3.769911184308
    effective_platelet_drag_nN_s_per_um: float = 85.2
    brownian_motion: bool = True
    remove_global_translation: bool = True
    adaptive_substeps: bool = True
    thermal_noise_model: str = "gaussian_independent"
    adaptive_target_displacement_um: float = 0.01
    max_internal_substeps: int = 10_000
    abort_displacement_um: float = 1.0

    # Rasterisation and mechanically generated spatial stack
    raster_resolution: int = 256
    clot_support_padding_um: float = 0.50
    mechanical_layers: int = 50
    mechanical_workers: int = 1
    layer_spacing_um: float = 0.070
    layer_seed_stride: int = 1
    save_intermediate_masks: bool = True

    # Void accessibility and optional first-passage experiment
    particle_radius_um: float = 0.10
    run_transport: bool = False
    tracer_diffusion_um2_s: float = 1.0
    transport_particles: int = 500
    transport_dt_s: float = 0.001
    transport_max_time_s: float = 5.0
    tracer_injection_rate_cells_s: float = 0.0

    def validate(self) -> None:
        positive = {
            "duration_s": self.duration_s,
            "dt_s": self.dt_s,
            "output_interval_s": self.output_interval_s,
            "domain_size_um": self.domain_size_um,
            "segment_target_um": self.segment_target_um,
            "temperature_K": self.temperature_K,
            "persistence_length_um": self.persistence_length_um,
            "contour_multiplier_C": self.contour_multiplier_C,
            "platelet_radius_um": self.platelet_radius_um,
            "platelet_force_radius_um": self.platelet_force_radius_um,
            "adaptive_target_displacement_um": self.adaptive_target_displacement_um,
            "raster_resolution": self.raster_resolution,
            "clot_support_padding_um": self.clot_support_padding_um,
            "mechanical_layers": self.mechanical_layers,
            "mechanical_workers": self.mechanical_workers,
            "layer_spacing_um": self.layer_spacing_um,
            "particle_radius_um": self.particle_radius_um,
        }
        for name, value in positive.items():
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be positive; received {value!r}")
        if self.equation_profile not in EQUATION_PROFILES:
            raise ValueError(f"equation_profile must be one of {EQUATION_PROFILES}")
        if self.parameter_profile not in PARAMETER_PROFILES:
            raise ValueError(f"parameter_profile must be one of {PARAMETER_PROFILES}")
        if self.numerical_profile not in NUMERICAL_PROFILES:
            raise ValueError(f"numerical_profile must be one of {NUMERICAL_PROFILES}")
        if self.geometry_profile not in GEOMETRY_PROFILES:
            raise ValueError(f"geometry_profile must be one of {GEOMETRY_PROFILES}")
        if self.mobility_model not in MOBILITY_MODELS:
            raise ValueError(f"mobility_model must be one of {MOBILITY_MODELS}")
        if self.thermal_noise_model not in THERMAL_NOISE_MODELS:
            raise ValueError(f"thermal_noise_model must be one of {THERMAL_NOISE_MODELS}")
        if self.filopodia_count_model not in FILOPODIA_COUNT_MODELS:
            raise ValueError(f"filopodia_count_model must be one of {FILOPODIA_COUNT_MODELS}")
        if self.filopodia_selection_rule not in FILOPODIA_SELECTION_RULES:
            raise ValueError(f"filopodia_selection_rule must be one of {FILOPODIA_SELECTION_RULES}")
        if self.filopodia_force_distribution not in FILOPODIA_FORCE_DISTRIBUTIONS:
            raise ValueError(
                f"filopodia_force_distribution must be one of {FILOPODIA_FORCE_DISTRIBUTIONS}"
            )
        if self.filopodia_handover_rule not in FILOPODIA_HANDOVER_RULES:
            raise ValueError(f"filopodia_handover_rule must be one of {FILOPODIA_HANDOVER_RULES}")
        if self.trajectory_detail not in TRAJECTORY_DETAILS:
            raise ValueError(f"trajectory_detail must be one of {TRAJECTORY_DETAILS}")
        if self.dt_s > self.duration_s:
            raise ValueError("dt_s cannot exceed duration_s")
        if not 0.0 < self.gamma_max < self.contour_multiplier_C:
            raise ValueError("gamma_max must satisfy 0 < gamma_max < C")
        if self.platelet_count < 0:
            raise ValueError("platelet_count cannot be negative")
        if not 0.0 <= self.platelet_near_node_fraction <= 1.0:
            raise ValueError("platelet_near_node_fraction must be in [0,1]")
        if not 0.0 < self.platelet_core_fraction <= 1.0:
            raise ValueError("platelet_core_fraction must be in (0,1]")
        if self.platelet_force_radius_um <= self.platelet_radius_um + self.fiber_radius_um:
            raise ValueError("platelet_force_radius_um must exceed the fibrin contact shell")
        if self.cohesion_search_stride < 1:
            raise ValueError("cohesion_search_stride must be at least 1")
        if self.cohesion_max_new_per_node_per_search < 0:
            raise ValueError("cohesion_max_new_per_node_per_search cannot be negative")
        if self.cohesion_max_total_neighbors < 0:
            raise ValueError("cohesion_max_total_neighbors cannot be negative")
        if self.filopodia_fixed_count < 0:
            raise ValueError("filopodia_fixed_count cannot be negative")
        if self.max_internal_substeps < 1:
            raise ValueError("max_internal_substeps must be at least 1")
        if self.layer_seed_stride < 1:
            raise ValueError("layer_seed_stride must be at least 1")
        if int(self.mechanical_workers) != self.mechanical_workers or self.mechanical_workers < 1:
            raise ValueError("mechanical_workers must be a positive integer")
        if self.run_transport:
            for name in ("tracer_diffusion_um2_s", "transport_dt_s", "transport_max_time_s"):
                if float(getattr(self, name)) <= 0.0:
                    raise ValueError(f"{name} must be positive")
            if self.transport_particles < 1:
                raise ValueError("transport_particles must be at least 1")
        if self.tracer_injection_rate_cells_s < 0.0:
            raise ValueError("tracer_injection_rate_cells_s cannot be negative")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, values: dict[str, Any]) -> "Program19Config":
        values = dict(values)
        # Program19 JSON files remain loadable. The old archive_literal label
        # is now the more precise codeocean_literal equation profile.
        if "physics_profile" in values and "equation_profile" not in values:
            old = values.pop("physics_profile")
            values["equation_profile"] = (
                "codeocean_literal" if old == "archive_literal" else old
            )
        known = {item.name for item in fields(cls)}
        unknown = sorted(set(values) - known)
        if unknown:
            raise ValueError(f"Unknown configuration keys: {', '.join(unknown)}")
        result = cls(**values)
        result.validate()
        return result

    @classmethod
    def load_json(cls, path: str | Path) -> "Program19Config":
        with Path(path).open("r", encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))

    def save_json(self, path: str | Path) -> None:
        with Path(path).open("w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2, ensure_ascii=False)
