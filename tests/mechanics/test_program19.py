from __future__ import annotations

import json
import tempfile
import threading
import unittest
from pathlib import Path

import numpy as np

from program19_core.bending import bending_triplet
from program19_core.cohesion import add_new_cohesion_links
from program19_core.config import Program19Config
from program19_core.experiments import paired_replicate_base_seed
from program19_core.integrator import mobility_and_diffusion
from program19_core.interactions import platelet_platelet_contact_forces
from program19_core.mechanical_runner import run_pipeline
from program19_core.continuation import continue_mechanical_run
from program19_core.catalog import scan_runs
from program19_core.evolving_transport import (
    StackHistory, load_stack_history, time_dependent_first_passage,
)
from program19_core.persistence import load_rng, load_state, save_checkpoint
from program19_core.stack_extension import extend_stack
from program19_core.profiles import make_profiled_config, profile_audit
from program19_core.state import SimulationState
from program19_core.transport import classify_void_components
from program19_core.through_clot import through_clot_first_passage
from program19_core.wlc import paper_stiffness_dF_dgamma_nN, wlc_force_nN


def empty_state(nodes, platelets) -> SimulationState:
    nodes = np.asarray(nodes, dtype=float).reshape((-1, 2))
    platelets = np.asarray(platelets, dtype=float).reshape((-1, 2))
    return SimulationState(
        nodes_um=nodes,
        platelets_um=platelets,
        backbone_pairs=np.empty((0, 2), dtype=int),
        backbone_L0_um=np.empty(0),
        bending_triples=np.empty((0, 3), dtype=int),
        bending_phi0_rad=np.empty(0),
        bending_dual_L0_um=np.empty(0),
        node_fiber_ids=[frozenset() for _ in nodes],
        incident_backbone=[[] for _ in nodes],
        original_adjacency=[[] for _ in nodes],
        filopodia=[[] for _ in platelets],
    )


class ForceLawTests(unittest.TestCase):
    def setUp(self):
        self.config = Program19Config(duration_s=0.001, dt_s=0.001, platelet_count=0)

    def test_wlc_reference_and_sign(self):
        L0 = 0.3
        self.assertAlmostEqual(wlc_force_nN(L0, L0, self.config), 0.0, places=12)
        self.assertGreater(wlc_force_nN(1.1 * L0, L0, self.config), 0.0)
        self.assertLess(wlc_force_nN(0.9 * L0, L0, self.config), 0.0)

    def test_wlc_tangent_matches_finite_difference(self):
        L0, gamma, epsilon = 0.3, 0.2, 1e-6
        numeric = (
            wlc_force_nN((1 + gamma + epsilon) * L0, L0, self.config)
            - wlc_force_nN((1 + gamma - epsilon) * L0, L0, self.config)
        ) / (2 * epsilon)
        analytic = paper_stiffness_dF_dgamma_nN(gamma, self.config)
        self.assertTrue(np.isclose(numeric, analytic, rtol=2e-5))

    def test_rest_angle_bending_zero_and_action_reaction(self):
        xi, xj, xk = np.asarray((0.0, 0.0)), np.asarray((1.0, 0.0)), np.asarray((2.0, 0.0))
        force = np.asarray(bending_triplet(xi, xj, xk, 0.0, 1.0, self.config))
        self.assertTrue(np.allclose(force, 0.0))
        bent = np.asarray(bending_triplet(xi, xj, np.asarray((2.0, 0.2)), 0.0, 1.0, self.config))
        self.assertTrue(np.allclose(np.sum(bent, axis=0), 0.0, atol=1e-12))

    def test_platelet_pair_is_reciprocal(self):
        state = empty_state([], [(0.0, 0.0), (2.1, 0.0)])
        force, overlaps, _ = platelet_platelet_contact_forces(state, self.config)
        self.assertEqual(overlaps, 1)
        self.assertTrue(np.allclose(force[0], -force[1]))


class ReproducibilityTests(unittest.TestCase):
    def test_corrected_free_stokes_diffusion(self):
        config = Program19Config(
            duration_s=0.001, dt_s=0.001, platelet_count=0, mobility_model="free_stokes"
        )
        _, _, node_D, platelet_D = mobility_and_diffusion(config)
        self.assertTrue(np.isclose(node_D, 1.098, rtol=0.01))
        self.assertTrue(np.isclose(platelet_D, 0.0486, rtol=0.01))

    def test_cohesion_zero_at_birth_and_valence(self):
        config = Program19Config(
            duration_s=0.001,
            dt_s=0.001,
            platelet_count=0,
            cohesion_distance_um=0.11,
            cohesion_max_total_neighbors=1,
            cohesion_max_new_per_node_per_search=1,
        )
        state = empty_state([(0.0, 0.0), (0.1, 0.0), (0.2, 0.0)], [])
        created = add_new_cohesion_links(state, config)
        self.assertEqual(created, 1)
        a, b = state.cohesion_pairs[0]
        self.assertAlmostEqual(
            state.cohesion_L0_um[0], np.linalg.norm(state.nodes_um[b] - state.nodes_um[a])
        )

    def test_mechanical_layer_controls_round_trip(self):
        config = Program19Config.from_dict(
            {
                "mechanical_layers": 4,
                "mechanical_workers": 3,
                "layer_seed_stride": 3,
                "save_intermediate_masks": True,
            }
        )
        self.assertEqual(config.mechanical_layers, 4)
        self.assertEqual(config.mechanical_workers, 3)
        self.assertEqual(config.layer_seed_stride, 3)
        self.assertTrue(config.save_intermediate_masks)

    def test_paired_replicates_use_non_overlapping_layer_seed_blocks(self):
        config = Program19Config(
            seed=2035,
            mechanical_layers=3,
            layer_seed_stride=1,
        )
        self.assertEqual(paired_replicate_base_seed(config, 0), 2035)
        self.assertEqual(paired_replicate_base_seed(config, 1), 2038)
        self.assertEqual(paired_replicate_base_seed(config, 2), 2041)

    def test_pore_classification_semantics(self):
        solid = np.ones((3, 5, 5), dtype=bool)
        solid[:, 1, 1] = False  # through-stack
        solid[0, 3, 1] = False  # blind z-connected
        solid[1, 3, 3] = False  # no z-face contact
        classes, summary = classify_void_components(solid, np.ones((5, 5), dtype=bool))
        self.assertEqual(classes[:, 1, 1].tolist(), [1, 1, 1])
        self.assertEqual(int(classes[0, 3, 1]), 2)
        self.assertEqual(int(classes[1, 3, 3]), 3)
        self.assertEqual(summary["component_count"], 3)

    def test_through_clot_transport_separates_transmission_return_and_retention(self):
        accessible = np.ones((5, 12, 12), dtype=bool)
        summary, times, outcomes, survival = through_clot_first_passage(
            accessible=accessible,
            dx_um=0.1,
            dz_um=0.1,
            diffusion_um2_s=0.5,
            particle_count=300,
            dt_s=0.001,
            max_time_s=0.5,
            rng=np.random.default_rng(17),
            injection_rates_cells_s=(10.0,),
            inlet_boundary="return",
        )
        self.assertEqual(len(times), 300)
        self.assertEqual(len(outcomes), 300)
        self.assertEqual(survival.shape, (101, 2))
        self.assertAlmostEqual(summary["probability_mass_check"], 1.0, places=14)
        self.assertGreater(summary["transmitted_count"], 0)
        self.assertGreater(summary["returned_count"], 0)
        scenario = summary["flux_scenarios"][0]
        self.assertAlmostEqual(
            scenario["transmission_rate_cells_s"],
            10.0 * summary["transmission_probability"],
        )

    def test_reflecting_inlet_has_no_returned_particles(self):
        accessible = np.ones((4, 8, 8), dtype=bool)
        summary, _, _, _ = through_clot_first_passage(
            accessible=accessible,
            dx_um=0.1,
            dz_um=0.1,
            diffusion_um2_s=1.0,
            particle_count=100,
            dt_s=0.001,
            max_time_s=0.5,
            rng=np.random.default_rng(3),
            inlet_boundary="reflecting",
        )
        self.assertEqual(summary["returned_count"], 0)
        self.assertAlmostEqual(summary["probability_mass_check"], 1.0, places=14)

    def test_smoke_run_has_exact_end_time_and_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = Program19Config(
                duration_s=0.0025,
                dt_s=0.001,
                output_interval_s=0.001,
                seed=7,
                domain_size_um=8.0,
                platelet_count=0,
                rho_degree3_um3=0.04,
                rho_degree4_um3=0.0,
                raster_resolution=32,
                mechanical_layers=2,
                mechanical_workers=2,
                layer_seed_stride=2,
                output_root=temporary,
                brownian_motion=False,
            )
            state, directory, status = run_pipeline(config)
            self.assertEqual(status, "completed")
            self.assertAlmostEqual(state.time_s, config.duration_s, places=14)
            manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["profiles"]["equation"], "paper_consistent")
            self.assertEqual(manifest["profiles"]["numerical"], "corrected_adaptive")
            self.assertEqual(manifest["execution"]["mechanical_workers_effective"], 2)
            self.assertEqual(manifest["stage_semantics"]["autoregression"], "not used")
            self.assertEqual(
                manifest["stage_semantics"]["interlayer_force_coupling"],
                "none; layers evolve independently",
            )
            self.assertTrue((directory / "pore_analysis.npz").exists())
            with np.load(directory / "mechanical_spatial_stack.npz") as stack:
                self.assertEqual(stack["solid_layers"].shape[0], 2)
                self.assertEqual(stack["layer_seeds"].tolist(), [7, 9])
                self.assertIn("no autoregression", str(stack["source"]))
            layer_dirs = sorted((directory / "mechanical_layers").glob("layer_*"))
            self.assertEqual(len(layer_dirs), 2)
            for layer_dir in layer_dirs:
                self.assertTrue((layer_dir / "checkpoint.npz").exists())
                with np.load(layer_dir / "trajectory.npz") as trajectory:
                    self.assertIn("cohesion_frame_offsets", trajectory.files)
                    self.assertIn("filopodia_node_ids", trajectory.files)
                with np.load(layer_dir / "intermediate_masks.npz") as masks:
                    self.assertEqual(masks["solid_masks"].shape[0], 4)

    def test_checkpoint_round_trip_preserves_state_and_rng(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = Program19Config(duration_s=0.01, dt_s=0.001, platelet_count=0)
            state = empty_state([(0.0, 0.0), (0.1, 0.0)], [])
            state.time_s = 0.004
            state.step = 4
            rng = np.random.default_rng(123)
            rng.normal(size=11)
            path = Path(temporary) / "checkpoint.npz"
            save_checkpoint(path, state, rng, config, source="test")
            expected = rng.normal(size=5)
            restored = load_state(path)
            restored_rng, policy = load_rng(path, fallback_seed=np.random.SeedSequence(9))
            self.assertEqual(policy, "exact_saved_rng_state")
            self.assertAlmostEqual(restored.time_s, state.time_s)
            self.assertTrue(np.allclose(restored.nodes_um, state.nodes_um))
            self.assertTrue(np.allclose(restored_rng.normal(size=5), expected))

    def test_parallel_layer_schedule_preserves_seeded_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            common = dict(
                duration_s=0.001,
                dt_s=0.001,
                output_interval_s=0.001,
                seed=71,
                domain_size_um=8.0,
                platelet_count=1,
                rho_degree3_um3=0.04,
                rho_degree4_um3=0.0,
                raster_resolution=32,
                mechanical_layers=2,
                run_transport=False,
            )
            serial = Program19Config(
                **common,
                mechanical_workers=1,
                output_root=str(Path(temporary) / "serial"),
            )
            parallel = Program19Config(
                **common,
                mechanical_workers=2,
                output_root=str(Path(temporary) / "parallel"),
            )
            _, serial_dir, _ = run_pipeline(serial)
            _, parallel_dir, _ = run_pipeline(parallel)
            with np.load(serial_dir / "mechanical_spatial_stack.npz") as saved:
                serial_stack = np.asarray(saved["solid_layers"])
                serial_seeds = np.asarray(saved["layer_seeds"])
            with np.load(parallel_dir / "mechanical_spatial_stack.npz") as saved:
                parallel_stack = np.asarray(saved["solid_layers"])
                parallel_seeds = np.asarray(saved["layer_seeds"])
            self.assertTrue(np.array_equal(serial_seeds, parallel_seeds))
            self.assertTrue(np.array_equal(serial_stack, parallel_stack))

    def test_parallel_continuation_preserves_seeded_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = Program19Config(
                duration_s=0.001,
                dt_s=0.001,
                output_interval_s=0.001,
                seed=91,
                domain_size_um=8.0,
                platelet_count=1,
                rho_degree3_um3=0.04,
                rho_degree4_um3=0.0,
                raster_resolution=24,
                mechanical_layers=2,
                mechanical_workers=1,
                output_root=str(Path(temporary) / "base"),
            )
            _, base, _ = run_pipeline(config)
            serial = continue_mechanical_run(
                base,
                additional_time_s=0.001,
                output_root=Path(temporary) / "continued_serial",
                mechanical_workers=1,
            )
            parallel = continue_mechanical_run(
                base,
                additional_time_s=0.001,
                output_root=Path(temporary) / "continued_parallel",
                mechanical_workers=2,
            )
            with np.load(serial / "mechanical_spatial_stack.npz") as saved:
                serial_stack = np.asarray(saved["solid_layers"])
                serial_seeds = np.asarray(saved["layer_seeds"])
            with np.load(parallel / "mechanical_spatial_stack.npz") as saved:
                parallel_stack = np.asarray(saved["solid_layers"])
                parallel_seeds = np.asarray(saved["layer_seeds"])
            self.assertTrue(np.array_equal(serial_seeds, parallel_seeds))
            self.assertTrue(np.array_equal(serial_stack, parallel_stack))

    def test_parallel_pipeline_can_start_from_gui_style_background_thread(self):
        with tempfile.TemporaryDirectory() as temporary:
            errors = []

            def target():
                try:
                    run_pipeline(
                        Program19Config(
                            duration_s=0.001,
                            dt_s=0.001,
                            output_interval_s=0.001,
                            seed=121,
                            domain_size_um=8.0,
                            platelet_count=0,
                            rho_degree3_um3=0.04,
                            rho_degree4_um3=0.0,
                            raster_resolution=24,
                            mechanical_layers=2,
                            mechanical_workers=2,
                            output_root=temporary,
                        )
                    )
                except Exception as exc:  # pragma: no cover - assertion reports detail
                    errors.append(exc)

            worker = threading.Thread(target=target)
            worker.start()
            worker.join(timeout=30.0)
            self.assertFalse(worker.is_alive())
            self.assertEqual(errors, [])

    def test_time_dependent_transport_matches_frozen_for_static_history(self):
        history = StackHistory(
            times_s=np.asarray([0.0, 0.05, 0.10]),
            solid_masks=np.zeros((3, 5, 12, 12), dtype=bool),
            analyzed_regions=np.ones((3, 5, 12, 12), dtype=bool),
            dx_um=0.1,
            dz_um=0.1,
            lineage=("synthetic",),
        )
        kwargs = dict(
            release_time_s=0.0,
            observation_time_s=0.1,
            particle_radius_um=0.01,
            diffusion_um2_s=0.5,
            particle_count=100,
            dt_s=0.001,
            recorded_particles=5,
        )
        dynamic, dynamic_arrays = time_dependent_first_passage(
            history, rng=np.random.default_rng(81), geometry_mode="evolving", **kwargs
        )
        frozen, frozen_arrays = time_dependent_first_passage(
            history, rng=np.random.default_rng(81), geometry_mode="frozen", **kwargs
        )
        self.assertEqual(dynamic["escaped_count"], frozen["escaped_count"])
        self.assertTrue(np.array_equal(dynamic_arrays["outcomes"], frozen_arrays["outcomes"]))
        self.assertTrue(np.allclose(dynamic_arrays["escape_times_s"], frozen_arrays["escape_times_s"], equal_nan=True))

    def test_frozen_transport_horizon_can_exceed_mechanical_history(self):
        history = StackHistory(
            times_s=np.asarray([0.0, 0.05, 0.10]),
            solid_masks=np.zeros((3, 5, 12, 12), dtype=bool),
            analyzed_regions=np.ones((3, 5, 12, 12), dtype=bool),
            dx_um=0.1,
            dz_um=0.1,
            lineage=("synthetic",),
        )
        summary, _ = time_dependent_first_passage(
            history,
            release_time_s=0.10,
            observation_time_s=0.30,
            particle_radius_um=0.01,
            diffusion_um2_s=0.5,
            particle_count=20,
            dt_s=0.001,
            rng=np.random.default_rng(82),
            geometry_mode="frozen",
            recorded_particles=0,
        )
        self.assertAlmostEqual(summary["observation_time_s"], 0.30)
        self.assertAlmostEqual(summary["mechanical_frame_time_s"], 0.10)

    def test_completed_run_can_continue_and_is_catalogued(self):
        with tempfile.TemporaryDirectory() as temporary:
            base_root = Path(temporary) / "base"
            continued_root = Path(temporary) / "continued"
            config = Program19Config(
                duration_s=0.001,
                dt_s=0.001,
                output_interval_s=0.001,
                seed=22,
                domain_size_um=8.0,
                platelet_count=0,
                rho_degree3_um3=0.04,
                rho_degree4_um3=0.0,
                raster_resolution=24,
                mechanical_layers=1,
                output_root=str(base_root),
                brownian_motion=False,
            )
            _, base, status = run_pipeline(config)
            self.assertEqual(status, "completed")
            continued = continue_mechanical_run(
                base, additional_time_s=0.001, output_root=continued_root
            )
            records = scan_runs((base_root, continued_root))
            self.assertEqual(len(records), 2)
            child = next(record for record in records if Path(record.path) == continued.resolve())
            self.assertAlmostEqual(child.mechanical_time_s, 0.002)
            self.assertEqual(Path(child.parent_path), base.resolve())

            extended = extend_stack(
                base,
                target_layer_count=2,
                first_new_seed=30,
                output_root=Path(temporary) / "extended",
            )
            with np.load(extended / "mechanical_spatial_stack.npz") as saved:
                self.assertEqual(saved["solid_layers"].shape[0], 2)
                self.assertEqual(saved["layer_seeds"].tolist(), [22, 30])
            history = load_stack_history(extended)
            self.assertEqual(history.solid_masks.shape[1], 2)

    def test_profiles_are_orthogonal_and_audited(self):
        config = make_profiled_config(
            equation_profile="paper_consistent",
            parameter_profile="codeocean_data_xml",
            numerical_profile="corrected_adaptive",
            geometry_profile="generated_2d",
        )
        self.assertEqual(config.equation_profile, "paper_consistent")
        self.assertEqual(config.platelet_fibrin_force_nN, 6.1)
        self.assertTrue(config.adaptive_substeps)
        self.assertFalse(config.cohesion_exclude_same_fiber)
        audit = profile_audit(config)
        self.assertEqual(audit["selected"]["parameter"], "codeocean_data_xml")

    def test_program19_json_alias_remains_loadable(self):
        config = Program19Config.from_dict({"physics_profile": "archive_literal"})
        self.assertEqual(config.equation_profile, "codeocean_literal")


if __name__ == "__main__":
    unittest.main()
