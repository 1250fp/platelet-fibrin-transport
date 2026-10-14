"""English Tk GUI for resumable Program19.3 mechanics and evolving transport."""

from __future__ import annotations

import json
import os
import queue
import threading
from pathlib import Path

from .config import (
    EQUATION_PROFILES,
    GEOMETRY_PROFILES,
    MOBILITY_MODELS,
    NUMERICAL_PROFILES,
    PARAMETER_PROFILES,
    Program19Config,
)
from .profiles import PROFILE_DESCRIPTIONS, make_profiled_config
from .mechanical_runner import run_pipeline


class Program19App:
    def __init__(self) -> None:
        import tkinter as tk
        from tkinter import ttk

        self.tk = tk
        self.ttk = ttk
        self.root = tk.Tk()
        self.root.title("Program19.3 — evolving clot mechanics and tracer transport")
        self.root.geometry("1120x850")
        self.messages: queue.Queue = queue.Queue()
        self.pause_event = threading.Event()
        self.pause_event.set()
        self.cancel_event = threading.Event()
        self.worker: threading.Thread | None = None
        self.last_run_directory: Path | None = None
        self.history_records = []
        self.variables: dict[str, object] = {}
        self._build()
        self.root.after(100, self._poll)

    def _var(self, key: str, value):
        if isinstance(value, bool):
            variable = self.tk.BooleanVar(value=value)
        else:
            variable = self.tk.StringVar(value=str(value))
        self.variables[key] = variable
        return variable

    def _entry(self, parent, row: int, label: str, key: str, default, width: int = 16):
        self.ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=8, pady=5)
        widget = self.ttk.Entry(parent, textvariable=self._var(key, default), width=width)
        widget.grid(row=row, column=1, sticky="ew", padx=8, pady=5)
        return widget

    def _build(self) -> None:
        defaults = make_profiled_config()
        outer = self.ttk.Frame(self.root, padding=10)
        outer.pack(fill="both", expand=True)
        notebook = self.ttk.Notebook(outer)
        notebook.pack(fill="both", expand=True)
        run = self.ttk.Frame(notebook, padding=12)
        profiles = self.ttk.Frame(notebook, padding=12)
        mechanics = self.ttk.Frame(notebook, padding=12)
        layers = self.ttk.Frame(notebook, padding=12)
        experiment = self.ttk.Frame(notebook, padding=12)
        history = self.ttk.Frame(notebook, padding=12)
        notebook.add(profiles, text="Profiles")
        notebook.add(run, text="Run")
        notebook.add(mechanics, text="Mechanics")
        notebook.add(layers, text="Mechanical layers and transport")
        notebook.add(experiment, text="Exploratory f0 sweep")
        notebook.add(history, text="Previous runs and transport")
        for tab in (profiles, run, mechanics, layers, experiment, history):
            tab.columnconfigure(1, weight=1)

        row = 0
        profile_specs = (
            ("Equation profile", "equation_profile", defaults.equation_profile, EQUATION_PROFILES),
            ("Parameter profile", "parameter_profile", defaults.parameter_profile, PARAMETER_PROFILES),
            ("Numerical profile", "numerical_profile", defaults.numerical_profile, NUMERICAL_PROFILES),
            ("Geometry profile", "geometry_profile", defaults.geometry_profile, GEOMETRY_PROFILES),
        )
        for label, key, value, choices in profile_specs:
            self.ttk.Label(profiles, text=label).grid(row=row, column=0, sticky="w", padx=8, pady=5)
            combo = self.ttk.Combobox(
                profiles,
                textvariable=self._var(key, value),
                values=choices,
                state="readonly",
            )
            combo.grid(row=row, column=1, sticky="ew", padx=8, pady=5)
            combo.bind("<<ComboboxSelected>>", lambda _event: self._update_profile_description())
            row += 1
        self.ttk.Button(
            profiles, text="Apply selected preset values", command=self._apply_profiles
        ).grid(row=row, column=0, columnspan=2, sticky="w", padx=8, pady=8); row += 1
        self.profile_description = self.tk.StringVar()
        self.ttk.Label(
            profiles,
            textvariable=self.profile_description,
            justify="left",
            wraplength=820,
        ).grid(row=row, column=0, columnspan=3, sticky="w", padx=8, pady=8)
        self._update_profile_description()

        row = 0
        self._entry(run, row, "Duration (s)", "duration_s", defaults.duration_s); row += 1
        self._entry(run, row, "Outer dt (s)", "dt_s", defaults.dt_s); row += 1
        self._entry(run, row, "Output interval (s)", "output_interval_s", defaults.output_interval_s); row += 1
        self._entry(run, row, "Random seed", "seed", defaults.seed); row += 1
        self._entry(run, row, "Platelet count", "platelet_count", defaults.platelet_count); row += 1
        self.ttk.Label(run, text="Mobility model").grid(row=row, column=0, sticky="w", padx=8, pady=5)
        self.ttk.Combobox(
            run,
            textvariable=self._var("mobility_model", defaults.mobility_model),
            values=MOBILITY_MODELS,
            state="readonly",
        ).grid(row=row, column=1, sticky="ew", padx=8, pady=5); row += 1
        self.ttk.Checkbutton(
            run, text="Brownian motion", variable=self._var("brownian_motion", defaults.brownian_motion)
        ).grid(row=row, column=0, sticky="w", padx=8, pady=5); row += 1
        self.ttk.Checkbutton(
            run,
            text="Remove global translation",
            variable=self._var("remove_global_translation", defaults.remove_global_translation),
        ).grid(row=row, column=0, sticky="w", padx=8, pady=5); row += 1
        self._entry(run, row, "Output root", "output_root", defaults.output_root, 48)
        self.ttk.Button(run, text="Browse…", command=self._browse_output).grid(row=row, column=2, padx=8); row += 1
        self._entry(run, row, "CodeOcean XML (optional)", "codeocean_xml_path", defaults.codeocean_xml_path, 48)
        self.ttk.Button(run, text="Browse…", command=self._browse_xml).grid(row=row, column=2, padx=8); row += 1

        row = 0
        self._entry(mechanics, row, "Contour multiplier C", "contour_multiplier_C", defaults.contour_multiplier_C); row += 1
        self._entry(mechanics, row, "Parallel chains N", "monomers_in_parallel_N", defaults.monomers_in_parallel_N); row += 1
        self._entry(mechanics, row, "Maximum strain γ", "gamma_max", defaults.gamma_max); row += 1
        self._entry(mechanics, row, "Persistence length (µm)", "persistence_length_um", defaults.persistence_length_um); row += 1
        self._entry(mechanics, row, "Bending EI (nN·µm²)", "bending_EI_nN_um2", defaults.bending_EI_nN_um2); row += 1
        self._entry(mechanics, row, "Cohesion search stride", "cohesion_search_stride", defaults.cohesion_search_stride); row += 1
        self._entry(mechanics, row, "Maximum cohesion valence", "cohesion_max_total_neighbors", defaults.cohesion_max_total_neighbors); row += 1
        self._entry(mechanics, row, "Filopodial baseline f₀ (nN)", "filopodia_f0_nN", defaults.filopodia_f0_nN); row += 1
        self._entry(mechanics, row, "Adaptive target move (µm)", "adaptive_target_displacement_um", defaults.adaptive_target_displacement_um); row += 1
        self._entry(mechanics, row, "Maximum internal substeps", "max_internal_substeps", defaults.max_internal_substeps); row += 1

        row = 0
        self.ttk.Label(
            layers,
            text=(
                "Each z layer is a separate explicit 2-D mechanical realization. "
                "No AR equation and no force coupling between layers are used."
            ),
            justify="left",
            wraplength=820,
        ).grid(row=row, column=0, columnspan=3, sticky="w", padx=8, pady=8); row += 1
        self._entry(layers, row, "Mechanical layer count N", "mechanical_layers", defaults.mechanical_layers); row += 1
        self._entry(
            layers, row, "Parallel layer processes (1 = serial)",
            "mechanical_workers", defaults.mechanical_workers,
        ); row += 1
        self.ttk.Label(
            layers,
            text=(
                f"Detected logical CPU threads: {os.cpu_count() or 1}. "
                "Each process evolves one independent layer; start with 2–4 and "
                "increase only if memory and cooling remain adequate."
            ),
            justify="left",
            wraplength=820,
        ).grid(row=row, column=0, columnspan=3, sticky="w", padx=8, pady=(0, 8)); row += 1
        self._entry(layers, row, "Layer seed stride", "layer_seed_stride", defaults.layer_seed_stride); row += 1
        self._entry(layers, row, "Raster resolution", "raster_resolution", defaults.raster_resolution); row += 1
        self._entry(layers, row, "Layer spacing (µm)", "layer_spacing_um", defaults.layer_spacing_um); row += 1
        self.ttk.Checkbutton(
            layers,
            text="Save intermediate positions and topology",
            variable=self._var("save_trajectory", defaults.save_trajectory),
        ).grid(row=row, column=0, sticky="w", padx=8, pady=5); row += 1
        self.ttk.Checkbutton(
            layers,
            text="Save intermediate raster masks",
            variable=self._var("save_intermediate_masks", defaults.save_intermediate_masks),
        ).grid(row=row, column=0, sticky="w", padx=8, pady=5); row += 1
        self._entry(layers, row, "Particle radius (µm)", "particle_radius_um", defaults.particle_radius_um); row += 1
        self.ttk.Checkbutton(
            layers,
            text="Run first-passage particles",
            variable=self._var("run_transport", defaults.run_transport),
        ).grid(row=row, column=0, sticky="w", padx=8, pady=5); row += 1
        self._entry(layers, row, "Transport particles", "transport_particles", defaults.transport_particles); row += 1
        self._entry(layers, row, "Tracer diffusion (µm²/s)", "tracer_diffusion_um2_s", defaults.tracer_diffusion_um2_s); row += 1
        self._entry(layers, row, "Transport maximum time (s)", "transport_max_time_s", defaults.transport_max_time_s); row += 1
        self._entry(
            layers,
            row,
            "Tracer injection J_in (cells/s; 0 leaves R_agg undefined)",
            "tracer_injection_rate_cells_s",
            defaults.tracer_injection_rate_cells_s,
        ); row += 1

        row = 0
        self.ttk.Checkbutton(
            experiment,
            text="Run an exploratory f0 sweep instead of one simulation",
            variable=self._var("experiment_enabled", False),
        ).grid(row=row, column=0, columnspan=2, sticky="w", padx=8, pady=8); row += 1
        self._entry(experiment, row, "f0 values (nN, comma separated)", "experiment_f0_values", "0,2.2,4.4,8.8", 38); row += 1
        self._entry(experiment, row, "Paired-seed replicates", "experiment_replicates", 1); row += 1
        self.ttk.Label(
            experiment,
            text=(
                "The sweep saves raw and aggregate CSV files. MET is the restricted mean "
                "first-passage time. R_agg is calculated only as J_in(1-P_escape); without "
                "a specified injection flux it remains undefined rather than being invented."
            ),
            justify="left",
            wraplength=820,
        ).grid(row=row, column=0, columnspan=3, sticky="w", padx=8, pady=12)

        self._build_history_tab(history)

        controls = self.ttk.Frame(outer)
        controls.pack(fill="x", pady=(10, 0))
        self.start_button = self.ttk.Button(controls, text="Start", command=self._start)
        self.start_button.pack(side="left", padx=4)
        self.pause_button = self.ttk.Button(controls, text="Pause", command=self._toggle_pause, state="disabled")
        self.pause_button.pack(side="left", padx=4)
        self.cancel_button = self.ttk.Button(controls, text="Cancel", command=self._cancel, state="disabled")
        self.cancel_button.pack(side="left", padx=4)
        self.ttk.Button(controls, text="Save config…", command=self._save_config).pack(side="left", padx=4)
        self.ttk.Button(controls, text="Load config…", command=self._load_config).pack(side="left", padx=4)
        self.open_button = self.ttk.Button(controls, text="Open results", command=self._open_results, state="disabled")
        self.open_button.pack(side="right", padx=4)
        self.progress = self.ttk.Progressbar(outer, maximum=100.0)
        self.progress.pack(fill="x", pady=(10, 4))
        self.status = self.tk.StringVar(value="Ready")
        self.ttk.Label(outer, textvariable=self.status).pack(anchor="w")
        self.root.after(250, self._refresh_history)

    def _build_history_tab(self, parent) -> None:
        legacy = Path(__file__).resolve().parents[1].parent / "program19_2" / "results"
        current = Path(__file__).resolve().parents[1] / "results"
        self.library_roots = self.tk.StringVar(value=f"{legacy};{current}")
        self.ttk.Label(parent, text="Result roots (semicolon separated)").grid(
            row=0, column=0, sticky="w", padx=6, pady=5
        )
        self.ttk.Entry(parent, textvariable=self.library_roots, width=90).grid(
            row=0, column=1, sticky="ew", padx=6, pady=5
        )
        self.ttk.Button(parent, text="Browse…", command=self._browse_library).grid(
            row=0, column=2, padx=6, pady=5
        )
        self.ttk.Button(parent, text="Refresh", command=self._refresh_history).grid(
            row=0, column=3, padx=6, pady=5
        )
        columns = ("created", "status", "time", "layers", "f0", "seed", "program")
        self.history_tree = self.ttk.Treeview(
            parent, columns=columns, show="headings", height=16, selectmode="browse"
        )
        headings = {
            "created": "Created",
            "status": "Status",
            "time": "Mechanical time (s)",
            "layers": "Layers",
            "f0": "f₀ (nN)",
            "seed": "Base seed",
            "program": "Program",
        }
        widths = {"created": 155, "status": 85, "time": 125, "layers": 75,
                  "f0": 80, "seed": 85, "program": 120}
        for key in columns:
            self.history_tree.heading(key, text=headings[key])
            self.history_tree.column(key, width=widths[key], anchor="center")
        self.history_tree.grid(row=1, column=0, columnspan=4, sticky="nsew", padx=6, pady=8)
        parent.rowconfigure(1, weight=1)
        scrollbar = self.ttk.Scrollbar(parent, orient="vertical", command=self.history_tree.yview)
        scrollbar.grid(row=1, column=4, sticky="ns", pady=8)
        self.history_tree.configure(yscrollcommand=scrollbar.set)

        buttons = self.ttk.Frame(parent)
        buttons.grid(row=2, column=0, columnspan=4, sticky="ew", padx=4, pady=4)
        self.ttk.Button(buttons, text="Inspect parameters", command=self._inspect_selected).pack(
            side="left", padx=4
        )
        self.ttk.Button(buttons, text="Continue mechanical evolution…", command=self._continue_selected).pack(
            side="left", padx=4
        )
        self.ttk.Button(buttons, text="Extend stack…", command=self._extend_selected).pack(
            side="left", padx=4
        )
        self.ttk.Button(buttons, text="Evolving transport…", command=self._evolving_selected).pack(
            side="left", padx=4
        )
        self.ttk.Button(buttons, text="Frozen transport…", command=self._frozen_selected).pack(
            side="left", padx=4
        )
        self.ttk.Button(buttons, text="Open selected folder", command=self._open_selected).pack(
            side="left", padx=4
        )
        self.ttk.Label(
            parent,
            text=(
                "Continuation preserves the saved mechanical state. Evolving transport uses "
                "the time-resolved masks; frozen transport reuses one selected saved frame for "
                "an observation horizon independent of the mechanical duration."
            ),
            justify="left",
            wraplength=980,
        ).grid(row=3, column=0, columnspan=4, sticky="w", padx=8, pady=8)

    def _browse_output(self) -> None:
        from tkinter import filedialog
        value = filedialog.askdirectory(initialdir=str(Path(__file__).resolve().parents[1]))
        if value:
            self.variables["output_root"].set(value)

    def _browse_library(self) -> None:
        from tkinter import filedialog
        value = filedialog.askdirectory(initialdir=str(Path(__file__).resolve().parents[1]))
        if value:
            existing = [item.strip() for item in self.library_roots.get().split(";") if item.strip()]
            if value not in existing:
                existing.append(value)
            self.library_roots.set(";".join(existing))
            self._refresh_history()

    def _refresh_history(self) -> None:
        from tkinter import messagebox
        from .catalog import scan_runs

        roots = [Path(item.strip()) for item in self.library_roots.get().split(";") if item.strip()]
        self.status.set("Scanning previous simulations…")
        self.root.update_idletasks()
        try:
            self.history_records = scan_runs(roots)
            for item in self.history_tree.get_children():
                self.history_tree.delete(item)
            for index, record in enumerate(self.history_records):
                self.history_tree.insert(
                    "", "end", iid=str(index), values=(
                        record.created_local.replace("T", " "),
                        record.status,
                        f"{record.mechanical_time_s:g}",
                        f"{record.completed_layers}/{record.requested_layers}",
                        f"{record.f0_nN:g}",
                        record.base_seed,
                        record.program,
                    ),
                )
            self.status.set(f"Found {len(self.history_records)} reusable mechanical runs")
        except Exception as exc:
            self.status.set("Could not scan previous simulations")
            messagebox.showerror("Simulation catalog", str(exc))

    def _selected_record(self):
        from tkinter import messagebox
        selected = self.history_tree.selection()
        if not selected:
            messagebox.showinfo("Previous simulations", "Select one simulation first.")
            return None
        try:
            return self.history_records[int(selected[0])]
        except (IndexError, ValueError):
            messagebox.showerror("Previous simulations", "The selected catalog row is no longer valid.")
            return None

    def _inspect_selected(self) -> None:
        from tkinter import messagebox
        record = self._selected_record()
        if record is None:
            return
        config_path = Path(record.path) / "config.json"
        config = json.loads(config_path.read_text(encoding="utf-8"))
        important = (
            "duration_s", "dt_s", "output_interval_s", "seed", "mechanical_layers",
            "mechanical_workers",
            "layer_spacing_um", "filopodia_f0_nN", "platelet_count", "mobility_model",
            "equation_profile", "parameter_profile", "numerical_profile", "geometry_profile",
            "particle_radius_um", "tracer_diffusion_um2_s", "transport_dt_s",
        )
        lines = [f"{key}: {config.get(key)}" for key in important if key in config]
        lines.extend((
            f"status: {record.status}",
            f"time history: {record.has_time_history}",
            f"exact RNG checkpoint: {record.exact_rng_checkpoint}",
            f"path: {record.path}",
        ))
        messagebox.showinfo("Simulation parameters", "\n".join(lines))

    def _prepare_background(self, text: str) -> bool:
        from tkinter import messagebox
        if self.worker is not None and self.worker.is_alive():
            messagebox.showinfo("Program19.3", "Another operation is still running.")
            return False
        self.cancel_event.clear()
        self.pause_event.set()
        self.start_button.configure(state="disabled")
        self.pause_button.configure(state="normal", text="Pause")
        self.cancel_button.configure(state="normal")
        self.open_button.configure(state="disabled")
        self.progress["value"] = 0
        self.status.set(text)
        return True

    def _continue_selected(self) -> None:
        from tkinter import filedialog, messagebox, simpledialog
        record = self._selected_record()
        if record is None:
            return
        if record.completed_layers != record.requested_layers:
            messagebox.showerror(
                "Mechanical continuation",
                "Temporal continuation currently requires every requested layer to be complete.",
            )
            return
        extra = simpledialog.askfloat(
            "Continue mechanical evolution",
            f"Current time: {record.mechanical_time_s:g} s\nAdditional physical time (s):",
            minvalue=1e-12,
            initialvalue=1.0,
            parent=self.root,
        )
        if extra is None:
            return
        workers = simpledialog.askinteger(
            "Continue mechanical evolution",
            "Maximum parallel layer processes (1 = serial):",
            minvalue=1,
            maxvalue=max(1, os.cpu_count() or 1),
            initialvalue=int(self.variables["mechanical_workers"].get()),
            parent=self.root,
        )
        if workers is None:
            return
        default_root = Path(__file__).resolve().parents[1] / "results"
        output = filedialog.askdirectory(
            title="Select Program19.3 output root", initialdir=str(default_root), mustexist=False
        )
        if not output or not self._prepare_background("Preparing temporal continuation…"):
            return

        def work() -> None:
            try:
                from .continuation import continue_mechanical_run
                directory = continue_mechanical_run(
                    record.path,
                    additional_time_s=extra,
                    output_root=output,
                    mechanical_workers=workers,
                    progress_callback=self.messages.put,
                    pause_event=self.pause_event,
                    cancel_event=self.cancel_event,
                )
                status_data = json.loads((Path(directory) / "status.json").read_text(encoding="utf-8"))
                self.messages.put({
                    "type": "finished", "status": status_data.get("status", "completed"),
                    "directory": str(directory),
                })
            except Exception as exc:
                self.messages.put({"type": "error", "error": str(exc)})

        self.worker = threading.Thread(target=work, daemon=True)
        self.worker.start()

    def _evolving_selected(self) -> None:
        from tkinter import filedialog, messagebox, simpledialog
        record = self._selected_record()
        if record is None:
            return
        if not record.has_time_history:
            messagebox.showerror(
                "Evolving transport",
                "This run does not contain intermediate_masks.npz for every layer.",
            )
            return
        config = Program19Config.load_json(Path(record.path) / "config.json")
        release = simpledialog.askfloat(
            "Evolving transport", "Tracer release time in the mechanical history (s):",
            initialvalue=0.0, minvalue=0.0, parent=self.root,
        )
        if release is None:
            return
        default_horizon = max(min(1.0, record.mechanical_time_s - release), config.transport_dt_s)
        horizon = simpledialog.askfloat(
            "Evolving transport",
            f"Observation time after release (s)\nAvailable through t={record.mechanical_time_s:g} s:",
            initialvalue=default_horizon, minvalue=config.transport_dt_s, parent=self.root,
        )
        if horizon is None:
            return
        radius = simpledialog.askfloat(
            "Evolving transport", "Tracer radius (µm):", initialvalue=config.particle_radius_um,
            minvalue=1e-12, parent=self.root,
        )
        if radius is None:
            return
        diffusion = simpledialog.askfloat(
            "Evolving transport", "Tracer diffusion coefficient (µm²/s):",
            initialvalue=config.tracer_diffusion_um2_s, minvalue=1e-12, parent=self.root,
        )
        if diffusion is None:
            return
        particles = simpledialog.askinteger(
            "Evolving transport", "Number of tracers:", initialvalue=config.transport_particles,
            minvalue=1, parent=self.root,
        )
        if particles is None:
            return
        include_frozen = messagebox.askyesno(
            "Paired comparison",
            "Also run the same Brownian tracers against the frozen release-time geometry?",
            parent=self.root,
        )
        default_root = Path(__file__).resolve().parents[1] / "results" / "transport_analyses"
        default_root.mkdir(parents=True, exist_ok=True)
        output = filedialog.askdirectory(
            title="Select transport-analysis output root", initialdir=str(default_root), mustexist=False
        )
        if not output or not self._prepare_background("Loading the time-resolved stack…"):
            return

        def work() -> None:
            try:
                from .evolving_transport import run_evolving_transport_comparison
                directory = run_evolving_transport_comparison(
                    record.path,
                    output_root=output,
                    release_time_s=release,
                    observation_time_s=horizon,
                    particle_radius_um=radius,
                    diffusion_um2_s=diffusion,
                    particle_count=particles,
                    dt_s=config.transport_dt_s,
                    seed=config.seed + 9303,
                    recorded_particles=min(25, particles),
                    trajectory_interval_s=max(config.transport_dt_s, 0.01),
                    include_frozen=include_frozen,
                    progress_callback=self.messages.put,
                )
                self.messages.put({"type": "finished", "status": "completed", "directory": str(directory)})
            except Exception as exc:
                self.messages.put({"type": "error", "error": str(exc)})

        self.worker = threading.Thread(target=work, daemon=True)
        self.worker.start()

    def _frozen_selected(self) -> None:
        from tkinter import filedialog, messagebox, simpledialog
        record = self._selected_record()
        if record is None:
            return
        if not record.has_time_history:
            messagebox.showerror(
                "Frozen transport",
                "This run does not contain intermediate_masks.npz for every layer.",
            )
            return
        config = Program19Config.load_json(Path(record.path) / "config.json")
        release = simpledialog.askfloat(
            "Frozen transport",
            f"Mechanical time of the frozen stack (s)\nAvailable: 0 to {record.mechanical_time_s:g} s:",
            initialvalue=record.mechanical_time_s,
            minvalue=0.0,
            maxvalue=record.mechanical_time_s,
            parent=self.root,
        )
        if release is None:
            return
        horizon = simpledialog.askfloat(
            "Frozen transport",
            "Observation time after release (s)\n"
            "This may exceed the saved mechanical duration:",
            initialvalue=config.transport_max_time_s,
            minvalue=config.transport_dt_s,
            parent=self.root,
        )
        if horizon is None:
            return
        radius = simpledialog.askfloat(
            "Frozen transport", "Tracer radius (µm):", initialvalue=config.particle_radius_um,
            minvalue=1e-12, parent=self.root,
        )
        if radius is None:
            return
        diffusion = simpledialog.askfloat(
            "Frozen transport", "Tracer diffusion coefficient (µm²/s):",
            initialvalue=config.tracer_diffusion_um2_s, minvalue=1e-12, parent=self.root,
        )
        if diffusion is None:
            return
        particles = simpledialog.askinteger(
            "Frozen transport", "Number of tracers:", initialvalue=config.transport_particles,
            minvalue=1, parent=self.root,
        )
        if particles is None:
            return
        default_root = Path(__file__).resolve().parents[1] / "results" / "transport_analyses"
        default_root.mkdir(parents=True, exist_ok=True)
        output = filedialog.askdirectory(
            title="Select frozen-transport output root",
            initialdir=str(default_root),
            mustexist=False,
        )
        if not output or not self._prepare_background("Loading the selected frozen stack…"):
            return

        def work() -> None:
            try:
                from .evolving_transport import run_frozen_transport_analysis
                directory = run_frozen_transport_analysis(
                    record.path,
                    output_root=output,
                    release_time_s=release,
                    observation_time_s=horizon,
                    particle_radius_um=radius,
                    diffusion_um2_s=diffusion,
                    particle_count=particles,
                    dt_s=config.transport_dt_s,
                    seed=config.seed + 9403,
                    recorded_particles=min(25, particles),
                    trajectory_interval_s=max(config.transport_dt_s, 0.01),
                    progress_callback=self.messages.put,
                )
                self.messages.put({
                    "type": "finished", "status": "completed", "directory": str(directory),
                })
            except Exception as exc:
                self.messages.put({"type": "error", "error": str(exc)})

        self.worker = threading.Thread(target=work, daemon=True)
        self.worker.start()

    def _extend_selected(self) -> None:
        from tkinter import filedialog, messagebox, simpledialog
        from .stack_extension import used_layer_seeds

        record = self._selected_record()
        if record is None:
            return
        if record.completed_layers != record.requested_layers:
            messagebox.showerror("Extend stack", "Only a complete stack can be extended.")
            return
        if record.parent_path and (Path(record.path) / "continuation_manifest.json").is_file():
            messagebox.showerror(
                "Extend stack",
                "Extend the original stack first; then continue the enlarged stack in time.",
            )
            return
        target = simpledialog.askinteger(
            "Extend stack",
            f"Current layer count: {record.completed_layers}\nNew total layer count:",
            initialvalue=record.completed_layers + 5,
            minvalue=record.completed_layers + 1,
            parent=self.root,
        )
        if target is None:
            return
        roots = [Path(item.strip()) for item in self.library_roots.get().split(";") if item.strip()]
        used = used_layer_seeds(roots)
        suggested = max(used, default=record.base_seed + record.completed_layers - 1) + 1
        first_seed = simpledialog.askinteger(
            "Extend stack",
            f"First seed for the {target - record.completed_layers} new layers:\n"
            f"Suggested collision-free value: {suggested}",
            initialvalue=suggested,
            minvalue=0,
            parent=self.root,
        )
        if first_seed is None:
            return
        workers = simpledialog.askinteger(
            "Extend stack",
            "Maximum parallel layer processes for the new layers (1 = serial):",
            minvalue=1,
            maxvalue=max(1, os.cpu_count() or 1),
            initialvalue=int(self.variables["mechanical_workers"].get()),
            parent=self.root,
        )
        if workers is None:
            return
        default_root = Path(__file__).resolve().parents[1] / "results"
        output = filedialog.askdirectory(
            title="Select Program19.3 output root", initialdir=str(default_root), mustexist=False
        )
        if not output or not self._prepare_background("Preparing stack extension…"):
            return

        def work() -> None:
            try:
                from .stack_extension import extend_stack
                directory = extend_stack(
                    record.path,
                    target_layer_count=target,
                    first_new_seed=first_seed,
                    output_root=output,
                    mechanical_workers=workers,
                    forbidden_seeds=used,
                    progress_callback=self.messages.put,
                )
                self.messages.put({"type": "finished", "status": "completed", "directory": str(directory)})
            except Exception as exc:
                self.messages.put({"type": "error", "error": str(exc)})

        self.worker = threading.Thread(target=work, daemon=True)
        self.worker.start()

    def _open_selected(self) -> None:
        record = self._selected_record()
        if record is not None and Path(record.path).exists():
            os.startfile(record.path)

    def _browse_xml(self) -> None:
        from tkinter import filedialog
        value = filedialog.askopenfilename(
            initialdir=str(Path(__file__).resolve().parents[1] / "reference_inputs"),
            filetypes=(("CodeOcean XML", "*.xml"), ("All files", "*.*")),
        )
        if value:
            self.variables["codeocean_xml_path"].set(value)

    def _update_profile_description(self) -> None:
        if not hasattr(self, "profile_description"):
            return
        selected = {
            "equation": self.variables["equation_profile"].get(),
            "parameter": self.variables["parameter_profile"].get(),
            "numerical": self.variables["numerical_profile"].get(),
            "geometry": self.variables["geometry_profile"].get(),
        }
        lines = []
        for axis, name in selected.items():
            item = PROFILE_DESCRIPTIONS[axis][name]
            lines.append(f"{axis.title()}: {item.title}. {item.purpose}")
        self.profile_description.set("\n\n".join(lines))

    def _apply_profiles(self) -> None:
        config = make_profiled_config(
            equation_profile=str(self.variables["equation_profile"].get()),
            parameter_profile=str(self.variables["parameter_profile"].get()),
            numerical_profile=str(self.variables["numerical_profile"].get()),
            geometry_profile=str(self.variables["geometry_profile"].get()),
        )
        values = config.to_dict()
        for key, variable in self.variables.items():
            if key in values:
                variable.set(values[key])
        self._update_profile_description()

    def _config(self) -> Program19Config:
        base = make_profiled_config(
            equation_profile=str(self.variables["equation_profile"].get()),
            parameter_profile=str(self.variables["parameter_profile"].get()),
            numerical_profile=str(self.variables["numerical_profile"].get()),
            geometry_profile=str(self.variables["geometry_profile"].get()),
        ).to_dict()
        integer_keys = {
            "seed", "platelet_count", "cohesion_search_stride",
            "cohesion_max_total_neighbors", "max_internal_substeps",
            "raster_resolution", "mechanical_layers", "mechanical_workers",
            "layer_seed_stride", "transport_particles",
        }
        boolean_keys = {
            "brownian_motion", "remove_global_translation", "save_trajectory",
            "save_intermediate_masks", "run_transport"
        }
        for key, variable in self.variables.items():
            value = variable.get()
            if key in integer_keys:
                base[key] = int(value)
            elif key in boolean_keys:
                base[key] = bool(value)
            elif key.startswith("experiment_"):
                continue
            elif key in (
                "equation_profile", "parameter_profile", "numerical_profile",
                "geometry_profile", "mobility_model", "output_root", "codeocean_xml_path",
            ):
                base[key] = str(value)
            else:
                base[key] = float(value)
        return Program19Config.from_dict(base)

    def _start(self) -> None:
        from tkinter import messagebox
        try:
            config = self._config()
        except Exception as exc:
            messagebox.showerror("Invalid configuration", str(exc))
            return
        self.cancel_event.clear()
        self.pause_event.set()
        self.start_button.configure(state="disabled")
        self.pause_button.configure(state="normal", text="Pause")
        self.cancel_button.configure(state="normal")
        self.open_button.configure(state="disabled")
        self.progress["value"] = 0
        self.status.set("Starting…")

        def work() -> None:
            try:
                if bool(self.variables["experiment_enabled"].get()):
                    from .experiments import run_f0_sweep

                    values = [
                        float(item.strip())
                        for item in str(self.variables["experiment_f0_values"].get()).split(",")
                        if item.strip()
                    ]
                    directory = run_f0_sweep(
                        config,
                        values,
                        int(self.variables["experiment_replicates"].get()),
                        progress_callback=self.messages.put,
                    )
                    status = "completed"
                else:
                    _, directory, status = run_pipeline(
                        config,
                        progress_callback=self.messages.put,
                        pause_event=self.pause_event,
                        cancel_event=self.cancel_event,
                    )
                self.messages.put({"type": "finished", "status": status, "directory": str(directory)})
            except Exception as exc:
                self.messages.put({"type": "error", "error": str(exc)})

        self.worker = threading.Thread(target=work, daemon=True)
        self.worker.start()

    def _toggle_pause(self) -> None:
        if self.pause_event.is_set():
            self.pause_event.clear()
            self.pause_button.configure(text="Resume")
            self.status.set("Paused")
        else:
            self.pause_event.set()
            self.pause_button.configure(text="Pause")

    def _cancel(self) -> None:
        self.cancel_event.set()
        self.pause_event.set()
        self.status.set("Cancelling safely…")

    def _poll(self) -> None:
        from tkinter import messagebox
        try:
            while True:
                payload = self.messages.get_nowait()
                kind = payload.get("type")
                if kind == "progress":
                    stage = payload.get("stage", "working")
                    fraction = float(payload.get("fraction", 0.0))
                    # Mechanics occupies 75% of the visual progress bar.
                    if stage in {"mechanics", "continuing_mechanics", "generating_added_layers"}:
                        value = 75.0 * fraction
                    elif stage in {"evolving_transport", "frozen_transport"}:
                        value = 100.0 * fraction
                    elif stage == "stack_assembly":
                        value = 90.0
                    elif stage == "pore_analysis":
                        value = 94.0
                    elif stage == "first_passage":
                        value = 97.0
                    elif stage in {
                        "mechanics", "continuing_mechanics", "generating_added_layers"
                    } and int(payload.get("parallel_workers", 1)) > 1:
                        self.status.set(
                            f"{stage.replace('_', ' ').title()} · "
                            f"{payload.get('completed_layers', 0)}/"
                            f"{payload.get('layer_count', 0)} layers completed · "
                            f"{payload.get('active_workers', 0)} processes active · "
                            f"latest L{payload.get('layer_index', 0)} "
                            f"t={payload.get('simulated_time_s', 0):.4g} s"
                        )
                    else:
                        value = 100.0
                    self.progress["value"] = value
                    if stage in {"evolving_transport", "frozen_transport"}:
                        self.status.set(
                            f"{stage.replace('_', ' ').title()} · "
                            f"t={payload.get('simulated_time_s', 0):.4g} s · "
                            f"active={payload.get('active_particles', 0)}"
                        )
                    else:
                        self.status.set(
                            f"{stage.replace('_', ' ').title()} · layer "
                            f"{payload.get('layer_index', payload.get('completed_layers', 0))}/"
                            f"{payload.get('layer_count', payload.get('requested_layers', 0))} · "
                            f"t={payload.get('simulated_time_s', 0):.4g} s "
                            f"· links={payload.get('cohesion_links', 0)}"
                        )
                elif kind == "finished":
                    self.last_run_directory = Path(payload["directory"])
                    self.status.set(f"{payload['status'].title()} · {self.last_run_directory}")
                    self.start_button.configure(state="normal")
                    self.pause_button.configure(state="disabled")
                    self.cancel_button.configure(state="disabled")
                    self.open_button.configure(state="normal")
                    self._refresh_history()
                    messagebox.showinfo("Program19.3", f"Run {payload['status']}\n\n{self.last_run_directory}")
                elif kind == "error":
                    self.status.set("Failed; partial results were preserved")
                    self.start_button.configure(state="normal")
                    self.pause_button.configure(state="disabled")
                    self.cancel_button.configure(state="disabled")
                    messagebox.showerror("Program19.3 failed", payload["error"])
        except queue.Empty:
            pass
        self.root.after(100, self._poll)

    def _save_config(self) -> None:
        from tkinter import filedialog, messagebox
        try:
            config = self._config()
            path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=(("JSON", "*.json"),))
            if path:
                config.save_json(path)
        except Exception as exc:
            messagebox.showerror("Cannot save configuration", str(exc))

    def _load_config(self) -> None:
        from tkinter import filedialog, messagebox
        path = filedialog.askopenfilename(filetypes=(("JSON", "*.json"),))
        if not path:
            return
        try:
            values = Program19Config.load_json(path).to_dict()
            for key, variable in self.variables.items():
                variable.set(values[key])
            self._update_profile_description()
        except Exception as exc:
            messagebox.showerror("Cannot load configuration", str(exc))

    def _open_results(self) -> None:
        if self.last_run_directory and self.last_run_directory.exists():
            os.startfile(self.last_run_directory)

    def run(self) -> None:
        self.root.mainloop()
