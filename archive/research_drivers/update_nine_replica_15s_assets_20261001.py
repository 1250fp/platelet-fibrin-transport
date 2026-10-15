"""Consolidate completed nine-replica frozen transport through mechanical 15 s.

Default: write only a new workspace analysis directory, never the manuscript.
--apply --backup PATH explicitly copies verified data/figures into the paper,
after checking a separately created, self-contained pre-edit backup.  No TEX
prose, manuscript PDF, bibliography, mechanics, transport, or ZIP is changed.
Do not run until the new batch has completed all 63 geometries / 252 chunks.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
from statistics import NormalDist

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/six_replica_manuscript_20260924"
NEW = ROOT / "outputs/transporte_nueve_replicas_15s_20261001"
OUT = ROOT / "outputs/nine_replica_15s_manuscript_20261001"
PAPER = Path('external_data/manuscript')
REPS = tuple(range(1, 10))
FORCES = (0.0, 2.2, 4.4, 8.8)
MECHANICAL_TIMES = (0, 5, 10, 15)
N = 256
HORIZON = 100.0
EXPECTED = dict(mechanical_replicas=9, nominal_states=144,
                unique_geometries=117, histories=29952,
                new_unique_geometries=63, new_chunks=252,
                added_unique_histories=16128, force_control_contrasts=81)
INPUTS: dict[str, str] = {}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def track(path):
    path = Path(path)
    require(path.is_file(), f"Missing required input: {path}")
    INPUTS[str(path)] = sha(path)
    return path


def load(path):
    return json.loads(track(path).read_text(encoding="utf-8-sig"))


def read_csv(path):
    with track(path).open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        for key, value in row.items():
            try:
                row[key] = float(value)
            except (ValueError, TypeError):
                pass
    return rows


def save(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n", encoding="utf-8")


def csvwrite(path, rows):
    require(bool(rows), f"No rows to write: {path}")
    keys = list(dict.fromkeys(k for row in rows for k in row))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def read_events(path):
    return np.atleast_1d(np.genfromtxt(track(path), delimiter=",", names=True))


def copy_events(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        f.write(",".join(data.dtype.names) + "\n")
        np.savetxt(f, np.column_stack([data[k] for k in data.dtype.names]),
                   delimiter=",", fmt="%.17g")


def wilson(k, n):
    z = NormalDist().inv_cdf(.975)
    p = k / n
    den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return center-half, center+half


def survival(data, times):
    times = np.asarray(times)
    return ((data["outcome"][:, None] == 0)
            | (data["time_s"][:, None] > times)).mean(axis=0)


def interval(delta, indices):
    return tuple(float(v) for v in
                 np.quantile(delta[indices].mean(axis=1), [.025, .975]))


def compare_existing(new_rows, old_rows, label):
    current = {r["state_id"]: r for r in new_rows}
    require(len(current) == len(new_rows), f"Duplicate {label} state identifiers")
    for old in old_rows:
        new = current[old["state_id"]]
        for key, value in old.items():
            require(key in new, f"Missing old field {label}/{key}")
            if isinstance(value, (float, int)):
                require(np.isclose(new[key], value, rtol=0, atol=1e-10),
                        f"Prior {label} changed: {old['state_id']} / {key}")
            else:
                require(new[key] == value,
                        f"Prior {label} changed: {old['state_id']} / {key}")


def check_protocol(old, new):
    # The new launcher records the same scalar names where possible; require
    # explicit protocol keys rather than guessing a changed scientific setup.
    wanted = {"radius_um": .25, "D_um2_s": .4, "dt_s": .000015625,
              "horizon_s": 100., "H_um": .1, "Lz_um": 5.}
    aliases = {"dt_s": ("transport_dt_s",), "horizon_s": ("T_s", "T"),
               "H_um": ("layer_spacing_um",), "Lz_um": ("thickness_um", "L_um")}
    for key, value in wanted.items():
        require(key in old and np.isclose(old[key], value), f"Old protocol {key}")
        found = next((new[k] for k in (key,) + aliases.get(key, ()) if k in new), None)
        require(found is not None and np.isclose(found, value),
                f"New protocol must explicitly preserve {key}={value}")


def validate_events(data, status, xyz, sid):
    require(len(data) == N and np.array_equal(data["particle"], np.arange(N)),
            f"Particle order/allocation invalid: {sid}")
    require(all(np.isfinite(data[k]).all() for k in data.dtype.names),
            f"Nonfinite transport output: {sid}")
    require(np.isin(data["outcome"], [0, 1, 2]).all(), f"Unknown outcome: {sid}")
    require(((data["time_s"] > 0) & (data["time_s"] <= HORIZON)).all(),
            f"Out-of-range event times: {sid}")
    require((data["time_s"][data["outcome"] == 0] == HORIZON).all(),
            f"Censoring must occur at {HORIZON}: {sid}")
    require((data["zf"][data["outcome"] == 1] == 0).all()
            and (data["zf"][data["outcome"] == 2] == 5).all(),
            f"Incorrect exit face: {sid}")
    require(len(status) == N and np.isin(status, [-1, 0, 1]).all(),
            f"Invalid connectivity classification: {sid}")
    require(np.array_equal(xyz, np.c_[data["x0"], data["y0"], data["z0"]]),
            f"Connectivity and transport starts differ: {sid}")
    require(not ((status == -1) & (data["outcome"] != 0)).any(),
            f"Certified-disconnected start escaped; stop, do not reinterpret: {sid}")


def event_metrics(data, status):
    escaped = data["outcome"] != 0
    rm = data["time_s"]
    lo, hi = wilson(int((~escaped).sum()), N)
    se = rm.std(ddof=1) / math.sqrt(N)
    nc, nd, nu = (int((status == k).sum()) for k in [1, -1, 0])
    return dict(escaped=int(escaped.sum()), censored=int((~escaped).sum()),
                S100=float((~escaped).mean()), S100_MC95_low=lo, S100_MC95_high=hi,
                RMST100_s=float(rm.mean()),
                RMST100_MC95_low=float(max(0., rm.mean()-1.96*se)),
                RMST100_MC95_high=float(min(HORIZON, rm.mean()+1.96*se)),
                F_conn_lower=nc/N, F_conn_upper=1-nd/N,
                connected=nc, disconnected=nd, unresolved=nu,
                escaped_from_independently_connected=int((escaped & (status == 1)).sum()),
                escaped_from_unresolved=int((escaped & (status == 0)).sum()),
                censored_from_connected=int((~escaped & (status == 1)).sum()),
                censored_from_disconnected=int((~escaped & (status == -1)).sum()),
                censored_from_unresolved=int((~escaped & (status == 0)).sum()),
                lateral_contact_particles=int((data["lateral_hits"] > 0).sum()),
                geometric_subdivisions=int(data["geometric_subdivisions"].sum()))


def paired_stats(a, b, indices):
    ds = 100 * ((a["outcome"] == 0).astype(float) - (b["outcome"] == 0).astype(float))
    dr = a["time_s"] - b["time_s"]
    sl, sh = interval(ds, indices)
    rl, rh = interval(dr, indices)
    xa = np.c_[a["x0"], a["y0"], a["z0"]]
    xb = np.c_[b["x0"], b["y0"], b["z0"]]
    return dict(delta_S100_pp=float(ds.mean()),
                delta_S100_MC95_low_pp=sl, delta_S100_MC95_high_pp=sh,
                delta_RMST100_s=float(dr.mean()),
                delta_RMST100_MC95_low_s=rl, delta_RMST100_MC95_high_s=rh,
                same_initial_at_same_particle_index=int((xa == xb).all(axis=1).sum()),
                shared_initial_coordinates_any_index=len(set(map(tuple, xa)) & set(map(tuple, xb))),
                delta_initial_z_mean_um=float(a["z0"].mean()-b["z0"].mean()))


def geometry_delta(row, base):
    return dict(delta_C_A_common_pp=100*(row["C_A_common"]-base["C_A_common"]),
                delta_C_A_layer_pp=100*(row["C_A_layer"]-base["C_A_layer"]),
                delta_accessible_fraction_pp=100*(row["phi_accessible_common"]-base["phi_accessible_common"]),
                accessible_volume_change_pct=100*(row["accessible_volume_um3"]/base["accessible_volume_um3"]-1))


def build_comparisons(rows, old_rows, out):
    metrics = ["S100", "RMST100_s", "C_A_common", "phi_accessible_common", "accessible_volume_um3"]
    aggregates = []
    for tm in MECHANICAL_TIMES:
        for force in ([0.] if tm == 0 else FORCES):
            selected = [r for r in rows if r["mechanical_time_s"] == tm and r["f0_nN"] == force]
            require(len(selected) == 9, f"Incomplete aggregate at {tm}/{force}")
            agg = dict(mechanical_time_s=tm, f0_nN=force, replicas=9)
            for metric in metrics:
                v = np.array([r[metric] for r in selected])
                agg.update({metric+"_mean": float(v.mean()), metric+"_min": float(v.min()),
                            metric+"_max": float(v.max()), metric+"_sd_between_replicas": float(v.std(ddof=1))})
            aggregates.append(agg)
    csvwrite(out/"aggregate_metrics.csv", aggregates)
    comparisons = []
    for tm in [0, 5, 10]:
        for force in ([0.] if tm == 0 else FORCES):
            before = [r for r in old_rows if r["mechanical_time_s"] == tm and r["f0_nN"] == force]
            after = [r for r in rows if r["mechanical_time_s"] == tm and r["f0_nN"] == force]
            c = dict(mechanical_time_s=tm, f0_nN=force, previous_replicas=6, current_replicas=9)
            for metric in metrics:
                a = float(np.mean([r[metric] for r in before]))
                b = float(np.mean([r[metric] for r in after]))
                c.update({metric+"_six": a, metric+"_nine": b, metric+"_cohort_change": b-a})
            comparisons.append(c)
    csvwrite(out/"comparison_six_vs_nine_common_times.csv", comparisons)
    temporal = []
    for force in FORCES:
        a = next(r for r in aggregates if r["mechanical_time_s"] == 10 and r["f0_nN"] == force)
        b = next(r for r in aggregates if r["mechanical_time_s"] == 15 and r["f0_nN"] == force)
        r = dict(f0_nN=force, replicas=9, time_before_s=10, time_after_s=15)
        for metric in metrics:
            r.update({metric+"_tm10": a[metric+"_mean"], metric+"_tm15": b[metric+"_mean"],
                      metric+"_temporal_change": b[metric+"_mean"]-a[metric+"_mean"]})
        temporal.append(r)
    csvwrite(out/"comparison_nine_tm10_vs_tm15.csv", temporal)
    lines = ["# Comparación numérica: seis y nueve réplicas, y extensión a 15 s", "",
             "Esta comparación separa dos cambios: ampliar la cohorte y avanzar el estado mecánico.", "",
             "## Cohorte en tiempos comunes", "",
             "| Tiempo mecánico (s) | f0 (nN) | S seis (%) | S nueve (%) | RMST seis (s) | RMST nueve (s) |",
             "|---:|---:|---:|---:|---:|---:|"]
    for r in comparisons:
        force = "Todas" if r["mechanical_time_s"] == 0 else f'{r["f0_nN"]:g}'
        lines.append(f'| {r["mechanical_time_s"]:g} | {force} | {100*r["S100_six"]:.2f} | '
                     f'{100*r["S100_nine"]:.2f} | {r["RMST100_s_six"]:.2f} | {r["RMST100_s_nine"]:.2f} |')
    lines += ["", "## Evolución mecánica con las mismas nueve réplicas", "",
              "| f0 (nN) | S tm10 (%) | S tm15 (%) | RMST tm10 (s) | RMST tm15 (s) |",
              "|---:|---:|---:|---:|---:|"]
    for r in temporal:
        lines.append(f'| {r["f0_nN"]:g} | {100*r["S100_tm10"]:.2f} | {100*r["S100_tm15"]:.2f} | '
                     f'{r["RMST100_s_tm10"]:.2f} | {r["RMST100_s_tm15"]:.2f} |')
    lines += ["", "Las medias dan el mismo peso a cada réplica mecánica. El transporte está congelado en cada estado; su horizonte sigue siendo 100 s.",
              "", "Los intervalos por contraste son bootstrap percentil de índices pareados, condicionales a cada red y al muestreo observado. No miden la incertidumbre entre coágulos ni el error numérico.",
              "", "Las nuevas réplicas no resuelven por sí mismas la equivalencia temporal inconclusa ni la conectividad no resuelta. Este archivo contiene comparaciones descriptivas, no una nueva declaración de validación.", ""]
    (out/"comparacion_numerica_seis_vs_nueve_y_15s.md").write_text("\n".join(lines), encoding="utf-8")
    return aggregates


def make_assets(rows, pairs, data, bykey, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figures = out/"figures"
    figures.mkdir(parents=True, exist_ok=True)
    lines = [r"\begin{table*}[t]", r"\centering\footnotesize",
             r"\setlength{\tabcolsep}{3.4pt}",
             r"\caption{Homogeneous frozen-state pilot. Values are means across nine paired mechanical replicas; bracketed transport ranges are minimum--maximum, not confidence intervals. Each unique geometry has 256 histories, tracer diameter $0.50~\mu\mathrm m$, $D=0.4~\mu\mathrm m^2/\mathrm s$, and $T_c=100$ s. The shared initial state is listed once. Contraction refers to the common analyzed region.}",
             r"\label{tab:homogeneous-pilot}", r"\begin{tabular}{@{}rrrrcc@{}}", r"\toprule",
             r"$t_m$ (s) & $f_0$ (nN) & $C_{A,\cap}$ (\%) & $\phi_{\rm acc}$ (\%) & $S(100)$ (\%) & RMST(100) (s) \\ \midrule"]
    for tm in MECHANICAL_TIMES:
        for force in ([0.] if tm == 0 else FORCES):
            rr = [bykey[rep, force, tm] for rep in REPS]
            ca = 100*np.mean([r["C_A_common"] for r in rr])
            phi = 100*np.mean([r["phi_accessible_common"] for r in rr])
            sv = 100*np.array([r["S100"] for r in rr])
            rm = np.array([r["RMST100_s"] for r in rr])
            force_text = "all" if tm == 0 else f"{force:g}"
            lines.append(f"{tm} & {force_text} & {ca:.3f} & {phi:.3f} & "
                         rf"${sv.mean():.2f}\ [{sv.min():.2f},{sv.max():.2f}]$ & "
                         rf"${rm.mean():.2f}\ [{rm.min():.2f},{rm.max():.2f}]$ " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    (out/"results_table.tex").write_text("\n".join(lines)+"\n", encoding="utf-8")
    plt.rcParams.update({"font.family": "cmr10", "mathtext.fontset": "cm", "font.size": 7.5,
                         "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
                         "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
    colors = dict(zip(REPS, ["#286e91", "#b95d2a", "#328765", "#865298", "#b48d16", "#555555",
                            "#ad3c62", "#447eab", "#526b39"]))
    markers = dict(zip(REPS, ["o", "s", "^", "D", "v", "P", "X", "<", "> ".strip()]))
    force_colors = {0.: "#4b6584", 2.2: "#207f91", 4.4: "#bf7838", 8.8: "#934c89"}

    def finish(fig, stem):
        for ext in ["pdf", "png"]:
            fig.savefig(figures/f"{stem}.{ext}", dpi=300)
        plt.close(fig)

    # Separate 10 and 15 s assets avoid silently changing an old figure's time.
    for tm in [10, 15]:
        fig, axs = plt.subplots(2, 1, figsize=(3.15, 4.6), layout="constrained", sharex=True)
        for ax, col, label, letter in zip(axs, ["C_A_common", "phi_accessible_common"],
                                        [r"$C_{A,\cap}$ (%)", r"$\phi_{\rm acc}$ (%)"], ["a", "b"]):
            for rep in REPS:
                ax.plot(FORCES, [100*bykey[rep, f, tm][col] for f in FORCES],
                        marker=markers[rep], color=colors[rep], ms=3.2, lw=.85, label=f"R{rep}")
            ax.set_ylabel(label)
            ax.set_xticks(FORCES)
            ax.text(.03, .94, letter, transform=ax.transAxes, weight="bold", fontsize=10, va="top")
        axs[-1].set_xlabel(r"$f_0$ (nN)")
        axs[-1].legend(frameon=False, loc="upper center", bbox_to_anchor=(.5, -.24),
                       fontsize=7, ncol=3, columnspacing=1.0, handlelength=1.5)
        finish(fig, f"homogeneous_geometry_tm{tm}")

        fig, axs = plt.subplots(2, 1, figsize=(3.15, 5.0), layout="constrained")
        curves = []
        for force, color in force_colors.items():
            dd = [data[bykey[rep, force, tm]["canonical_geometry"]] for rep in REPS]
            tt = np.unique(np.r_[0., HORIZON, *[d["time_s"] for d in dd]])
            yy = np.array([survival(d, tt) for d in dd])
            axs[0].fill_between(tt, yy.min(axis=0), yy.max(axis=0), step="post", color=color, alpha=.09, lw=0)
            axs[0].step(tt, yy.mean(axis=0), where="post", color=color, lw=1.1, label=f"{force:g} nN")
            for j, rep in enumerate(REPS):
                curves.extend(dict(f0_nN=force, replicate=rep, mechanical_time_s=tm,
                                   t_s=float(t), S=float(s)) for t, s in zip(tt, yy[j]))
        csvwrite(out/f"survival_curves_tm{tm}.csv", curves)
        axs[0].set(xlabel="Transport time (s)", ylabel=r"$S(t)$", xlim=(0, 100), ylim=(0, 1.04))
        axs[0].legend(frameon=False, fontsize=7, ncol=2, loc="lower left")
        for rep in REPS:
            rr = [p for p in pairs if p["replicate"] == rep and p["mechanical_time_s"] == tm]
            xx = np.array([p["f0_nN"] for p in rr]) + (rep-5)*.10
            yy = np.array([p["delta_RMST100_s"] for p in rr])
            lo = np.array([p["delta_RMST100_MC95_low_s"] for p in rr])
            hi = np.array([p["delta_RMST100_MC95_high_s"] for p in rr])
            axs[1].errorbar(xx, yy, yerr=[yy-lo, hi-yy], fmt=markers[rep]+"-", color=colors[rep],
                            ms=3.2, lw=.75, capsize=1.6, label=f"R{rep}")
        axs[1].axhline(0, color=".4", ls="--", lw=.7)
        axs[1].set(xlabel=r"$f_0$ (nN)", ylabel=r"$\Delta$RMST(100 s) (s)", xticks=FORCES[1:])
        axs[1].legend(frameon=False, fontsize=7, ncol=3, loc="upper center",
                       bbox_to_anchor=(.5, -.24), columnspacing=1., handlelength=1.5)
        for ax, letter in zip(axs, ["a", "b"]):
            ax.text(.95, .94, letter, transform=ax.transAxes, weight="bold", fontsize=10, va="top")
        finish(fig, f"homogeneous_transport_tm{tm}")
    # Keep legacy filenames at their original mechanical time; TEX must opt in
    # explicitly to the additional 15 s assets.
    for stem in ["homogeneous_geometry", "homogeneous_transport"]:
        for ext in ["pdf", "png"]:
            shutil.copy2(figures/f"{stem}_tm10.{ext}", figures/f"{stem}.{ext}")
    fig, axs = plt.subplots(2, 1, figsize=(3.15, 4.55), layout="constrained", sharex=True)
    for ax, metric, factor, label, letter in zip(axs, ["S100", "RMST100_s"], [100, 1],
                                               [r"$S(100)$ (%)", "RMST(100 s) (s)"], ["a", "b"]):
        for force, color in force_colors.items():
            yy = factor*np.array([[bykey[rep, force, tm][metric] for tm in MECHANICAL_TIMES] for rep in REPS])
            ax.fill_between(MECHANICAL_TIMES, yy.min(axis=0), yy.max(axis=0), color=color, alpha=.09, lw=0)
            ax.plot(MECHANICAL_TIMES, yy.mean(axis=0), "o-", color=color, ms=3.4, lw=1, label=f"{force:g} nN")
        ax.set_ylabel(label)
        ax.text(.03, .94, letter, transform=ax.transAxes, weight="bold", fontsize=10, va="top")
    axs[0].legend(frameon=False, fontsize=7, ncol=2, loc="best")
    axs[-1].set(xlabel="Mechanical time (s)", xticks=MECHANICAL_TIMES)
    finish(fig, "homogeneous_time_evolution")
    save(out/"figure_manifest.json", dict(
        legacy_filenames_mechanical_time_s=10, new_explicit_endpoint_s=15,
        homogeneous_geometry="Replica curves; symbols identify nine independent seed groups.",
        homogeneous_transport="Mean S(t) with min--max range; paired per-replica RMST intervals are conditional 95% bootstrap intervals.",
        homogeneous_time_evolution="Frozen-state means and min--max ranges; not moving-boundary transport.",
        visual_QA_required=True))


def apply_assets(out, backup):
    backup = Path(backup).resolve()
    require(backup != PAPER.resolve(), "Backup must be separate from the active paper")
    for name in ["main2vJPS.tex", "main2vJPS.pdf", "homogeneous_pilot_sections.tex",
                 "validation_supplement.tex", "validation_supplement.pdf", "refs.bib",
                 "figures/refs.bib", "data/homogeneous_pilot/state_metrics.csv"]:
        require((backup/name).is_file(), f"Incomplete pre-edit backup: {backup/name}")
    for name in ["main2vJPS.tex", "main2vJPS.pdf", "homogeneous_pilot_sections.tex",
                 "validation_supplement.tex", "validation_supplement.pdf"]:
        active = PAPER/name
        if name == "main2vJPS.pdf" and not active.exists():
            active = PAPER/"6replica.pdf"
        require(sha(backup/name) == sha(active),
                f"Active paper differs from pre-edit backup: {name}; apply before TEX edits")
    require(int(load(backup/"data/homogeneous_pilot/resumen.json")["mechanical_replicas"]) == 6,
            "The comparison backup must contain the preceding six-replica paper")
    zip_path = PAPER/"minipaper_overleaf.zip"
    zip_before = sha(zip_path) if zip_path.exists() else None
    for p in sorted(out.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(out)
        dest = (PAPER/rel if rel.parts[0] == "figures"
                else PAPER/"data/homogeneous_pilot"/rel)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
    save(PAPER/"data/homogeneous_pilot/provenance.json", dict(
        consolidated_analysis=str(out), manifest_sha256=sha(out/"manifest.json"),
        previous_manuscript_snapshot=str(backup),
        new_transport_executed_in_source_batch=True, simulations_executed_by_this_script=False,
        tex_and_manuscript_pdfs_not_updated_by_this_script=True))
    require(zip_before is None or sha(zip_path) == zip_before, "ZIP changed unexpectedly")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="Explicitly copy completed data/figures into E manuscript")
    parser.add_argument("--backup", type=Path, help="Verified autocontained six-replica manuscript backup; required by --apply")
    args = parser.parse_args()
    require(not args.apply or args.backup is not None, "--apply requires --backup")
    execution = load(NEW/"execution_summary.json")
    require(execution.get("completed") is True, "New transport is not complete; no partial analysis permitted")
    jobs = execution["jobs"]
    require(len(jobs) == EXPECTED["new_chunks"], "Expected exactly 252 completed new chunks")
    old_protocol, new_protocol = load(BASE/"protocol.json"), load(NEW/"protocol.json")
    check_protocol(old_protocol, new_protocol)
    for path, digest in new_protocol["source_sha256"].items():
        require(sha(track(path)) == digest, f"Protected source changed: {path}")
    for job in jobs:
        require(job["exit_code"] == 0, f"Transport chunk failed: {job['job']['name']}")
        require(sha(track(job["job"]["output"])) == job["sha256"], "Chunk hash mismatch")
        require(job["job"]["N"] == 64 and job["job"]["dt"] == .000015625
                and job["job"]["T"] == HORIZON, "Unexpected chunk transport protocol")
    old_rows = read_csv(BASE/"state_metrics.csv")
    old_pairs = read_csv(BASE/"contrastes_vs_control.csv")
    old_geometry = read_csv(BASE/"common_geometry.csv")
    new_states = load(NEW/"states.json")
    require(len(old_rows) == 72 and len(new_states) == 72, "Expected 72 prior and 72 new labels")
    old_ids = {r["state_id"] for r in old_rows}
    require(not old_ids.intersection(s["state_id"] for s in new_states), "New states duplicate finished prior labels")
    new_unique = {s["canonical_geometry"] for s in new_states}
    require(len(new_unique) == 63, "Expected 63 new canonical geometries")
    require(set(j["job"]["state_id"] for j in jobs) == new_unique, "Chunk states do not match new geometry inventory")
    all_states = old_rows + new_states
    expected_keys = {(rep, f, tm) for rep in REPS for f in FORCES for tm in MECHANICAL_TIMES}
    actual_keys = {(int(r["replicate"]), float(r["f0_nN"]), int(r["mechanical_time_s"])) for r in all_states}
    require(len(all_states) == 144 and actual_keys == expected_keys, "Cohort key set incomplete or duplicated")
    canonical = sorted({r["canonical_geometry"] for r in all_states})
    require(len(canonical) == 117, "Canonical count mismatch")
    for r in all_states:
        if r["mechanical_time_s"] == 0:
            require(r["canonical_geometry"] == f"rep{int(r['replicate'])}_f0_0_tm0", "Incorrect t0 reuse")
        else:
            require(r["canonical_geometry"] == r["state_id"], "Postinitial states must be unique")

    # Baseline area is reconstructed from the original archive, never reset to
    # a continuation's first frame. Check every new state's stated contraction.
    initial_common, initial_mean = {}, {}
    for rep in REPS:
        if rep <= 6:
            rr = [r for r in old_geometry if r["replicate"] == rep and r["mechanical_time_s"] == 0]
            initial_common[rep] = float(rr[0]["common_area_um2"])
        else:
            s = next(s for s in new_states if s["replicate"] == rep and s["mechanical_time_s"] == 0)
            initial_common[rep] = float(s["common_area_um2"])
            initial_mean[rep] = float(s["mean_support_area_um2"])
    # Read old state inventories to recover the original mean layer support.
    old_batches = [ROOT/"outputs/analisis_homogeneo_dominio_comun_20260919",
                   ROOT/"outputs/transporte_replicas_3_4_20260923",
                   ROOT/"outputs/transporte_replica5_20260924", ROOT/"outputs/transporte_replica6_20260924"]
    for batch in old_batches:
        p = batch/"states.json"
        if p.exists():
            for s in load(p):
                if s["mechanical_time_s"] == 0 and "mean_support_area_um2" in s:
                    initial_mean[int(s["replicate"])] = float(s["mean_support_area_um2"])
    for rep in (1, 2):
        source = ROOT/"outputs/analisis_homogeneo_con_cohesion_20260917/states"/f"rep{rep}_f0_0_tm0"/"raster_with_cohesion.npz"
        with np.load(track(source), allow_pickle=False) as a:
            region, dx = a["analyzed_region_masks"], float(a["dx_um"])
        initial_mean[rep] = float(region.sum(axis=(1, 2)).mean()*dx**2)
    require(set(initial_mean) == set(REPS), "Missing original t0 mean layer areas; do not normalize to continuation start")
    geometry = list(old_geometry)
    for s in new_states:
        rep = int(s["replicate"])
        require(np.isclose(s["C_A_common"], 1-s["common_area_um2"]/initial_common[rep], atol=1e-10, rtol=0),
                f"Common contraction is not original-t0 normalized: {s['state_id']}")
        require(np.isclose(s["C_A_layer"], 1-s["mean_support_area_um2"]/initial_mean[rep], atol=1e-10, rtol=0),
                f"Layer contraction is not original-t0 normalized: {s['state_id']}")
        with np.load(track(NEW/"states"/s["canonical_geometry"]/"volume_samples.npz")) as v:
            inside = v["inside_common"]
            n = int(inside.sum())
            nv = int((inside & v["void"]).sum())
            na = int((inside & v["accessible"]).sum())
            require(n > 0 and len(inside) == 131072, "Wrong geometric MC proposal allocation")
        al, ah = wilson(na, n)
        vl, vh = wilson(nv, n)
        geometry.append(dict(state_id=s["state_id"], replicate=rep, f0_nN=s["f0_nN"],
                             mechanical_time_s=s["mechanical_time_s"], C_A_layer=s["C_A_layer"],
                             C_A_common=s["C_A_common"], common_area_um2=s["common_area_um2"],
                             common_volume_um3=5*s["common_area_um2"], phi_void_common=nv/n,
                             phi_accessible_common=na/n, phi_acc_MC95_low=al, phi_acc_MC95_high=ah,
                             phi_void_MC95_low=vl, phi_void_MC95_high=vh,
                             accessible_volume_um3=5*s["common_area_um2"]*na/n,
                             MC_proposals=131072, MC_points_in_domain=n, MC_accessible_points=na,
                             cohesive_segments=s["cohesive_segments"]))
    geom = {r["state_id"]: r for r in geometry}
    data, metrics, status_by_id = {}, {}, {}
    connrows, audits, disconnected = [], [], []
    for sid in canonical:
        if sid in new_unique:
            source = NEW/"states"/sid
            selected = sorted((j for j in jobs if j["job"]["state_id"] == sid), key=lambda j: j["job"]["offset"])
            require([j["job"]["offset"] for j in selected] == [0, 64, 128, 192], f"Chunk offsets incomplete: {sid}")
            d = np.concatenate([read_events(j["job"]["output"]) for j in selected])
            with np.load(track(source/"connectivity/initial_connectivity.npz")) as c:
                status, xyz = c["status"].copy(), c["xyz_um"].copy()
            conn = load(source/"connectivity/summary.json")
        else:
            source = BASE/"states"/sid
            d = read_events(source/"transport_merged.csv")
            cc = read_csv(source/"connectivity_escape.csv")
            require([r["particle"] for r in cc] == list(range(N)), f"Old connectivity particle order: {sid}")
            status = np.array([r["connectivity_before_transport"] for r in cc], dtype=int)
            xyz = np.c_[d["x0"], d["y0"], d["z0"]]
            require(all(r["escape_event"] == d["outcome"][i] and r["observed_time_s"] == d["time_s"][i]
                        for i, r in enumerate(cc)), f"Old event/classification archive differs: {sid}")
            conn = {}
        validate_events(d, status, xyz, sid)
        m = event_metrics(d, status)
        for k in ["connected", "disconnected", "unresolved", "F_conn_lower", "F_conn_upper"]:
            require(k not in conn or np.isclose(conn[k], m[k]), f"Connectivity count/bound mismatch {sid}/{k}")
        metrics[sid], data[sid], status_by_id[sid] = m, d, status
        connrows.append(dict(state_id=sid, **{k: m[k] for k in
                        ["connected", "disconnected", "unresolved", "F_conn_lower", "F_conn_upper"]}))
        audits.append(dict(state_id=sid, n=N, censored=m["censored"],
                           lateral_contacts=int(d["lateral_hits"].sum()), particles_touching_lateral=m["lateral_contact_particles"]))
        for i in np.flatnonzero(status == -1):
            disconnected.append(dict(state_id=sid, particle_zero_based=int(i), xyz_um=xyz[i].tolist(),
                                     outcome=int(d["outcome"][i]), time_s=float(d["time_s"][i])))
    rows = []
    for s in all_states:
        g = geom[s["state_id"]]
        rows.append(dict(state_id=s["state_id"], canonical_geometry=s["canonical_geometry"],
                         replicate=int(s["replicate"]), f0_nN=float(s["f0_nN"]), mechanical_time_s=float(s["mechanical_time_s"]),
                         **{k: g[k] for k in ["C_A_layer", "C_A_common", "phi_accessible_common", "accessible_volume_um3"]},
                         **metrics[s["canonical_geometry"]]))
    compare_existing(rows, old_rows, "state metrics")
    old_byid = {r["state_id"]: r for r in old_rows}
    rows = [dict(old_byid[r["state_id"]]) if r["state_id"] in old_byid else r for r in rows]
    rows.sort(key=lambda r: (r["replicate"], r["mechanical_time_s"], r["f0_nN"]))
    bykey = {(int(r["replicate"]), float(r["f0_nN"]), int(r["mechanical_time_s"])): r for r in rows}
    rng = np.random.default_rng(20260919)
    indices = {rep: rng.integers(0, N, size=(10000, N)) for rep in REPS}
    pairs, horizons, changes, times, temporal_pairs = [], [], [], [], []
    for rep in REPS:
        for tm in [5, 10, 15]:
            control = bykey[rep, 0., tm]
            b = data[control["canonical_geometry"]]
            for force in FORCES[1:]:
                r = bykey[rep, force, tm]
                a = data[r["canonical_geometry"]]
                pairs.append(dict(state_id=r["state_id"], control_state_id=control["state_id"], replicate=rep,
                                  mechanical_time_s=tm, f0_nN=force, **geometry_delta(r, control),
                                  **paired_stats(a, b, indices[rep])))
                for t in [10, 25, 50, 75, 100]:
                    horizons.append(dict(state_id=r["state_id"], replicate=rep, f0_nN=force, mechanical_time_s=tm,
                                         horizon_s=t, delta_S_pp=float(100*(survival(a, [t])[0]-survival(b, [t])[0])),
                                         delta_RMST_s=float(np.minimum(a["time_s"], t).mean()-np.minimum(b["time_s"], t).mean())))
        for force in FORCES:
            base = bykey[rep, force, 0]
            for tm in [5, 10, 15]:
                r = bykey[rep, force, tm]
                changes.append(dict(state_id=r["state_id"], replicate=rep, f0_nN=force, mechanical_time_s=tm,
                                    delta_C_A_common_pp=100*(r["C_A_common"]-base["C_A_common"]),
                                    delta_accessible_fraction_pp=100*(r["phi_accessible_common"]-base["phi_accessible_common"]),
                                    delta_S100_pp=100*(r["S100"]-base["S100"]), delta_RMST100_s=r["RMST100_s"]-base["RMST100_s"]))
            r, prev = bykey[rep, force, 15], bykey[rep, force, 10]
            temporal_pairs.append(dict(state_id=r["state_id"], reference_state_id=prev["state_id"],
                                       replicate=rep, f0_nN=force, from_time_s=10, to_time_s=15,
                                       **geometry_delta(r, prev),
                                       **paired_stats(data[r["canonical_geometry"]], data[prev["canonical_geometry"]], indices[rep])))
    require(len(pairs) == 81, "Contrast count mismatch")
    compare_existing(pairs, old_pairs, "paired bootstrap contrasts")
    old_pair_byid = {r["state_id"]: r for r in old_pairs}
    pairs = [dict(old_pair_byid[r["state_id"]]) if r["state_id"] in old_pair_byid else r for r in pairs]
    for r in rows:
        d = data[r["canonical_geometry"]]
        for t in [10, 25, 50, 75, 100]:
            times.append(dict(state_id=r["state_id"], T_s=t, S=float(survival(d, [t])[0]),
                              RMST_s=float(np.minimum(d["time_s"], t).mean())))

    # All input/result consistency checks have passed before any analysis output.
    OUT.mkdir(parents=True, exist_ok=True)
    for sid in canonical:
        d, status = data[sid], status_by_id[sid]
        target = OUT/"states"/sid
        copy_events(target/"transport_merged.csv", d)
        csvwrite(target/"connectivity_escape.csv", [dict(particle=i, connectivity_before_transport=int(status[i]),
                 escape_event=int(d["outcome"][i]), observed_time_s=float(d["time_s"][i])) for i in range(N)])
    for name, values in [("state_metrics.csv", rows), ("common_geometry.csv", geometry),
                         ("contrastes_vs_control.csv", pairs), ("cambios_desde_inicio.csv", changes),
                         ("sensibilidad_al_horizonte.csv", horizons), ("survival_horizons.csv", times),
                         ("auditoria_trayectorias.csv", audits), ("connectivity_summary.csv", connrows),
                         ("contrastes_tm15_vs_tm10.csv", temporal_pairs)]:
        csvwrite(OUT/name, values)
    csvwrite(OUT/"cohort_inventory.csv", [dict(state_id=r["state_id"], canonical_geometry=r["canonical_geometry"],
             replicate=r["replicate"], f0_nN=r["f0_nN"], mechanical_time_s=r["mechanical_time_s"],
             source_cohort=str(BASE if r["state_id"] in old_ids else NEW)) for r in rows])
    sums = {k: sum(m[k] for m in metrics.values()) for k in
            ["connected", "disconnected", "unresolved", "censored", "censored_from_connected",
             "censored_from_disconnected", "censored_from_unresolved", "lateral_contact_particles"]}
    require(sums["connected"]+sums["disconnected"]+sums["unresolved"] == 29952, "Classification total mismatch")
    summary = dict(**EXPECTED, **sums, total_force_control_contrasts=len(pairs),
                   positive_geometry_contrasts=sum(p["delta_C_A_common_pp"] > 0 and p["delta_accessible_fraction_pp"] > 0 for p in pairs),
                   min_accessible_volume_change_pct=min(p["accessible_volume_change_pct"] for p in pairs),
                   max_accessible_volume_change_pct=max(p["accessible_volume_change_pct"] for p in pairs),
                   maximum_same_index_initial_count=max(p["same_initial_at_same_particle_index"] for p in pairs),
                   minimum_shared_initial_count=min(p["shared_initial_coordinates_any_index"] for p in pairs),
                   maximum_shared_initial_count=max(p["shared_initial_coordinates_any_index"] for p in pairs),
                   total_lateral_events=sum(a["lateral_contacts"] for a in audits),
                   certified_disconnected=disconnected)
    summary["by_mechanical_time"] = {}
    for tm in MECHANICAL_TIMES:
        unique_ids = {r["canonical_geometry"] for r in rows if r["mechanical_time_s"] == tm}
        pp = [p for p in pairs if p["mechanical_time_s"] == tm]
        summary["by_mechanical_time"][str(tm)] = dict(unique_geometries=len(unique_ids), histories=N*len(unique_ids),
            **{k: sum(metrics[sid][k] for sid in unique_ids) for k in sums},
            active_control_contrasts=len(pp),
            negative_RMST_point_contrasts=sum(p["delta_RMST100_s"] < 0 for p in pp),
            RMST_conditional_intervals_below_zero=sum(p["delta_RMST100_MC95_high_s"] < 0 for p in pp),
            positive_contraction_and_accessibility=sum(p["delta_C_A_common_pp"] > 0 and p["delta_accessible_fraction_pp"] > 0 for p in pp))
    save(OUT/"resumen.json", summary)
    save(OUT/"consistency_checks.json", dict(passed=True,
         previous_six_replica_state_metrics_and_bootstrap_intervals_unchanged=True,
         new_chunk_hashes_verified=len(jobs), protected_source_files_verified=len(new_protocol["source_sha256"]),
         contraction_normalized_to_original_t0_for_all_new_states=True, disconnected_but_escaped=[],
         new_unique_geometries=63, unique_initial_states_counted_once=True))
    for source, name in [(BASE, "previous_six_replica_consolidation"), (NEW, "new_replicas_and_15s")]:
        for filename in ["protocol.json", "execution_summary.json"]:
            save(OUT/"source_protocols"/name/filename, load(source/filename))
    shutil.copytree(BASE/"source_protocols", OUT/"source_protocols"/"historical_batches", dirs_exist_ok=True)
    save(OUT/"protocol.json", dict(status="completed; reuse of prior six-replica data plus completed new batch",
         cohorts=[str(BASE), str(NEW)], mechanical_replicas=9, mechanical_times_s=list(MECHANICAL_TIMES),
         N_per_unique_geometry=N, nominal_states=144, unique_geometries=117, unique_histories=29952,
         radius_um=.25, D_um2_s=.4, dt_s=.000015625, horizon_s=100., H_um=.1, Lz_um=5.,
         all_replicas_prespecified_before_first_pilot=False, bootstrap_resamples=10000,
         bootstrap_seed=20260919, bootstrap_replica_draw_order=list(REPS),
         cohort_extension_sequential=True, simulations_executed_by_this_script=False,
         original_t0_contraction_reference_preserved=True,
         uncertainty="Marginal paired-index intervals conditional on each network; network ranges are not CIs."))
    save(OUT/"execution_summary.json", dict(completed=True, consolidation_only=True,
         prior_unique_histories=13824, added_unique_histories=16128, total_unique_histories=29952,
         simulations_executed_by_this_script=False, new_transport_source=str(NEW)))
    build_comparisons(rows, old_rows, OUT)
    make_assets(rows, pairs, data, bykey, OUT)
    save(OUT/"provenance.json", dict(previous_six_replica_analysis=str(BASE), new_transport_source=str(NEW),
         simulations_executed_by_this_script=False, manuscript_edits_not_performed_by_this_script=True))
    save(OUT/"manifest.json", dict(source_sha256=INPUTS, script_sha256=sha(__file__),
         analysis="Nine-replica completed-data consolidation at mechanical 0/5/10/15 s",
         output_sha256={str(p.relative_to(OUT)): sha(p) for p in sorted(OUT.rglob("*"))
                        if p.is_file() and p.name != "manifest.json"},
         bootstrap_resamples=10000, bootstrap_seed=20260919, replica_draw_order=list(REPS),
         visual_QA_required=True, temporal_equivalence_not_reassessed=True))
    if args.apply:
        apply_assets(OUT, args.backup)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("Assets applied to manuscript" if args.apply else "Workspace-only outputs; manuscript untouched")


if __name__ == "__main__":
    main()
