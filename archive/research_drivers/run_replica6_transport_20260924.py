"""Extend the frozen reflecting pilot to mechanical replica 6.

No mechanical integration, manuscript changes, ZIP update, or outcome analysis.
The old executables and protocol are reused; preparation is audited against one
previously analyzed state. Completed chunks are reusable; partial chunks fail.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ['MPLBACKEND'] = 'Agg'
sys.dont_write_bytecode = True
import numpy as np
import common_connectivity_analysis as conn
import run_reflection_100s as transport

ROOT = Path('C:/Users/1250fp/Documents/Codex/2026-07-14/es')
EXP = Path('E:/tesisprog/fibrina/sub_fibrin_model/program19_3/results/exploratory_experiments')
OLD = ROOT / 'outputs/piloto_homogeneo_reflexion_20260919'
GEOM = ROOT / 'outputs/analisis_homogeneo_con_cohesion_20260917'
OUT = ROOT / 'outputs/transporte_replica6_20260924'
BATCHES = [EXP / 'replica06_mechanics_20260924_035846']
REPLICAS = (6,)
NOMINAL_STATES, UNIQUE_STATES = 12*len(REPLICAS), 9*len(REPLICAS)
TRANSPORT_JOBS = 4*UNIQUE_STATES
TIMES = (0., 5., 10.)
H, L, RADIUS, DIFFUSION, DT, HORIZON = .1, 5., .25, .4, .000015625, 100.
N, POOL, MC = 256, 65536, 131072
protected = {}


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def protect(p):
    p = Path(p)
    protected[str(p)] = sha(p)
    return p


def load(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def save(p, data):
    p = Path(p)
    temp = p.with_suffix(p.suffix + '.tmp')
    temp.write_text(json.dumps(data, indent=2, allow_nan=False), encoding='utf-8')
    temp.replace(p)


def csvsave(p, rows):
    with Path(p).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        w.writeheader()
        w.writerows(rows)


def status(stage, **kw):
    save(OUT / 'status.json', dict(status=stage, pid=os.getpid(), updated_local=time.strftime('%Y-%m-%dT%H:%M:%S'), **kw))


def candidates(rep):
    rng = np.random.default_rng([916002, rep])
    angle = rng.uniform(0, 2*np.pi, POOL)
    rho = 1.19*np.sqrt(rng.random(POOL))
    z = rng.uniform(2, 3, POOL)
    return np.c_[rho*np.cos(angle), rho*np.sin(angle), z]


def query(scene, points, output, cap=.4):
    subprocess.run([str(OUT/'source_snapshot/common_clearance.exe'), 'points', str(scene), str(points), str(output), str(cap)],
                   check=True, capture_output=True)
    d = np.loadtxt(output)
    assert np.isfinite(d).all() and ((d >= 0) & (d <= cap)).all()
    return d


def segment_rows(nodes, pairs, radius, k):
    a, b = nodes[pairs[:, 0]], nodes[pairs[:, 1]]
    return np.c_[a, b, np.full(len(a), radius), np.full(len(a), k*H), np.full(len(a), (k+1)*H)]


def rebuild(layer_paths, times, raster):
    buffers = {t: [] for t in times}
    links = {t: [] for t in times}
    regions = {t: [] for t in times}
    solids = {t: [] for t in times}
    metadata = {t: [] for t in times}
    assert len(layer_paths) == 50
    for k, path in enumerate(layer_paths):
        cfg = load(protect(path/'config.json'))
        assert cfg['raster_resolution'] == 256 and cfg['domain_size_um'] == 47.6
        assert cfg['fiber_radius_um'] == .05
        summary = load(protect(path/'layer_summary.json'))
        assert summary['status'] == 'completed' and summary['error'] is None
        assert abs(summary['simulated_time_s']-10) < 1e-8
        with np.load(protect(path/'trajectory.npz'), allow_pickle=False) as tr:
            ts = tr['times_s']
            assert abs(ts[-1]-10) < 1e-8 and abs(ts[0]) < 1e-12
            nodes, platelets, backbone = tr['nodes_um'], tr['platelets_um'], tr['backbone_pairs']
            offsets, cpflat = tr['cohesion_frame_offsets'], tr['cohesion_pairs_flat']
            births = tr['cohesion_birth_times_s_flat']
            with np.load(protect(path/'intermediate_masks.npz'), allow_pickle=False) as masks:
                oldsolid, oldregion = masks['solid_masks'], masks['analyzed_region_masks']
            for t in times:
                j = int(np.argmin(np.abs(ts-t)))
                assert abs(ts[j]-t) < 1e-8
                cp = cpflat[offsets[j]:offsets[j+1]].reshape(-1, 2)
                assert np.all(births[offsets[j]:offsets[j+1]] <= t+1e-8)
                pairs = {tuple(sorted(map(int, pair))) for pair in cp}
                assert len(pairs) == len(cp)
                assert not pairs.intersection(tuple(sorted(map(int, pair))) for pair in backbone)
                nd, pl = nodes[j], platelets[j]
                fibre = segment_rows(nd, backbone, cfg['fiber_radius_um'], k)
                disks = np.c_[pl, pl, np.full(len(pl), cfg['platelet_radius_um']), np.full(len(pl), k*H), np.full(len(pl), (k+1)*H)]
                buffers[t].append(np.vstack([fibre, disks]))
                links[t].append(segment_rows(nd, cp, cfg['fiber_radius_um'], k))
                solid, region, extent, dx = raster(nd, pl, np.vstack([backbone, cp]), SimpleNamespace(**cfg))
                assert not np.any(oldsolid[j] & ~solid)
                assert not np.any(oldregion[j] & ~region)
                if len(cp) == 0:
                    assert np.array_equal(solid, oldsolid[j]) and np.array_equal(region, oldregion[j])
                regions[t].append(region)
                solids[t].append(solid)
                metadata[t].append(dict(layer=k+1, seed=cfg['seed'], path=str(path), frame_index=j,
                    saved_time_s=float(ts[j]), backbone_primitives=len(backbone), platelets=len(pl),
                    cohesive_links_saved=len(cp), cohesion_is_geometric_obstacle=True, cohesive_radius_um=cfg['fiber_radius_um']))
    return {t: dict(scene=np.vstack([np.concatenate(buffers[t]), np.concatenate(links[t])]),
                    region=np.stack(regions[t]), solid=np.stack(solids[t]), metadata=metadata[t], extent=extent, dx=dx)
            for t in times}


def write_geometry(dest, g):
    dest.mkdir(parents=True, exist_ok=False)
    np.savetxt(dest/'solid_full.tsv', g['scene'], fmt='%.17g')
    common = g['region'].all(axis=0)
    assert common.any()
    np.savez_compressed(dest/'geometry.npz', analyzed_region_masks=g['region'], solid_layers=g['solid'],
                        common_support=common, extent_um=g['extent'], dx_um=g['dx'], H_um=H, Lz_um=L)
    conn.write_domain(dest/'geometry.npz', dest/'domain.txt')
    return common


def prepare():
    assert not (OUT/'protocol.json').exists(), 'Preparation already complete; use run'
    OUT.mkdir(exist_ok=True)
    assert not (OUT/'source_snapshot').exists(), 'Partial preparation exists; inspect before retry'
    status('preparing', phase='protocol_and_source_audit')
    old = load(OLD/'protocol.json')
    assert load(OLD/'execution_summary.json')['completed']
    assert all(sha(p) == h for p, h in old['source_sha256'].items()), 'Prior pilot sources changed'
    assert old['N'] == N and old['T'] == HORIZON and old['dt_s'] == DT
    assert old['radius_um'] == RADIUS and old['D_um2_s'] == DIFFUSION and old['H_um'] == H
    snap = OUT/'source_snapshot'
    snap.mkdir()
    for name in [Path(__file__).name, 'run_reflection_100s.py', 'common_connectivity_analysis.py',
                 'transport_reflection_long.cpp', 'transport_reflection_v2.cpp', 'transport_geometry_bridge.cpp',
                 'common_domain_connectivity.cpp', 'homogeneous_clearance.cpp']:
        shutil.copy2(ROOT/'work'/name, snap/name)
    binaries = {'transport_reflection_long.exe': ROOT/'outputs/verificacion_reflexion_100s_20260919/transport_reflection_long.exe',
                'common_connectivity.exe': OLD/'common_connectivity.exe',
                'common_clearance.exe': ROOT/'outputs/analisis_homogeneo_dominio_comun_20260919/common_clearance.exe'}
    geom_manifest = load(ROOT/'outputs/analisis_homogeneo_dominio_comun_20260919/manifest.json')
    for name, source in binaries.items():
        hashes = old['source_sha256'] if name != 'common_clearance.exe' else geom_manifest['source_sha256']
        assert sha(source) == hashes[str(source)]
        shutil.copy2(source, snap/name)
    core = BATCHES[0]/'source_snapshot/program19_core'
    shutil.copytree(core, snap/'program19_core', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for p in snap.rglob('*'):
        if p.is_file(): protect(p)
    sys.path.insert(0, str(snap))
    raster = importlib.import_module('program19_core.layers').rasterize_geometry
    prior_prep = load(GEOM/'preparation_manifest.json')
    for name in ('layers.py', 'config.py', 'state.py'):
        oldpath = 'E:\\tesisprog\\fibrina\\sub_fibrin_model\\program19_3\\program19_core\\'+name
        assert sha(core/name) == prior_prep['source_sha256'][oldpath]
    for folder in ('states', 'initials', 'shared', 'transport', 'regression'):
        (OUT/folder).mkdir()

    # Replay an established geometry, including the exact primitive ordering.
    sid = 'rep1_f0_2.2_tm10'
    prior = load(GEOM/'states'/sid/'state.json')
    g = rebuild([Path(x['path']) for x in prior['layers_source']], [10.], raster)[10.]
    check = OUT/'regression'/sid
    write_geometry(check, g)
    assert sha(check/'solid_full.tsv') == sha(GEOM/'states'/sid/'solid_full.tsv')
    assert sha(check/'domain.txt') == sha(OLD/'states'/sid/'domain.txt')
    p = check/'candidates.tsv'
    np.savetxt(p, candidates(1), fmt='%.17g')
    distance = query(check/'solid_full.tsv', p, check/'clearance.tsv')
    chosen = candidates(1)[np.flatnonzero(distance >= RADIUS+1e-12)[:N]]
    assert np.array_equal(chosen, np.loadtxt(GEOM/'states'/sid/'initial_batch00.tsv'))
    save(OUT/'regression/preparation_replay.json', dict(passed=True, reference_state=sid,
        continuous_scene_byte_identical=True, common_domain_byte_identical=True, initial_coordinates_identical=True))
    print('REGRESSION: scene, domain and 256 initial coordinates match the previous pilot.', flush=True)

    allstates, inventory, unique, seedsets = [], [], {}, {}
    baseline_areas = {}
    for rep in REPLICAS:
        np.savetxt(OUT/'shared'/f'release_candidates_rep{rep}.tsv', candidates(rep), fmt='%.17g')
        rng = np.random.default_rng(np.random.SeedSequence([919201, rep]))
        xyz = rng.random((MC, 3))*[47.6, 47.6, 5]+[-23.8, -23.8, 0]
        np.savetxt(OUT/'shared'/f'volume_rep{rep}.tsv', xyz, fmt='%.17g')
    for batch in BATCHES:
        progress = load(protect(batch/'progress.json'))
        assert progress['status'] == 'completed' and progress['completed_stacks'] == 4
        plan = load(protect(batch/'replica_manifest.json'))
        rep = plan['replicate_id']
        assert rep in REPLICAS and plan['run_transport'] is False
        seedsets[rep] = plan['layer_seeds']
        assert len(seedsets[rep]) == len(set(seedsets[rep])) == 50
        for rel, digest in plan['source_snapshot_sha256'].items():
            assert sha(batch/'source_snapshot'/rel) == digest
        for name in ('layers.py', 'config.py', 'state.py'):
            assert sha(batch/'source_snapshot/program19_core'/name) == sha(core/name)
        runs = load(protect(batch/'completed_runs.json'))
        assert sorted(q['f0_nN'] for q in runs) == [0, 2.2, 4.4, 8.8]
        for row in runs:
            run = Path(row['run_directory'])
            cfg = load(protect(run/'config.json'))
            assert row['status'] == 'completed' and not row['transport_run']
            assert cfg['duration_s'] == 10 and cfg['output_interval_s'] == .1 and cfg['layer_spacing_um'] == H
            assert cfg['filopodia_f0_nN'] == row['f0_nN'] and cfg['mechanical_layers'] == 50
            paths = [run/'mechanical_layers'/f'layer_{k+1:03d}_seed{seed}' for k, seed in enumerate(seedsets[rep])]
            status('preparing', phase='geometry', replicate=rep, f0_nN=row['f0_nN'], completed_states=len(allstates))
            stack = rebuild(paths, TIMES, raster)
            for g in stack.values():
                assert [m['seed'] for m in g['metadata']] == seedsets[rep]
            for t, g in stack.items():
                sid = f"rep{rep}_f0_{row['f0_nN']:g}_tm{t:g}"
                dest = OUT/'states'/sid
                common = write_geometry(dest, g)
                scenehash, domhash = sha(dest/'solid_full.tsv'), sha(dest/'domain.txt')
                key = (rep, scenehash, domhash)
                canonical = unique.setdefault(key, sid)
                area = float(common.sum()*g['dx']**2)
                meanarea = float(g['region'].sum(axis=(1, 2)).mean()*g['dx']**2)
                if rep not in baseline_areas:
                    assert t == 0
                    baseline_areas[rep] = (area, meanarea)
                st = dict(state_id=sid, replicate=rep, f0_nN=row['f0_nN'], mechanical_time_s=t,
                    canonical_geometry=canonical, unique_geometry=(canonical == sid), source_run=str(run),
                    scene=str(dest/'solid_full.tsv'), domain=str(dest/'domain.txt'), initial=str(dest/'initial.tsv'),
                    seed=919010+100*(rep-1), common_area_um2=area, mean_support_area_um2=meanarea,
                    C_A_common=1-area/baseline_areas[rep][0], C_A_layer=1-meanarea/baseline_areas[rep][1],
                    layers_source=g['metadata'], solid_sha256=scenehash, domain_sha256=domhash,
                    cohesive_segments=sum(x['cohesive_links_saved'] for x in g['metadata']))
                allstates.append(st)
                save(dest/'state.json', st)
                inventory.append({k: v for k, v in st.items() if k != 'layers_source'})
            print(f'GEOMETRY: replica {rep}, f0={row["f0_nN"]:g}; 0/5/10 s rebuilt.', flush=True)
    used_prior_seeds = set(range(4000, 4100))
    for prior_batch in (EXP/'replica03_mechanics_20260922_151219', EXP/'replica04_mechanics_20260922_232615', EXP/'replica05_mechanics_20260923_173648'):
        used_prior_seeds.update(load(protect(prior_batch/'replica_manifest.json'))['layer_seeds'])
    assert set(seedsets) == set(REPLICAS)
    for rep in REPLICAS:
        assert set(seedsets[rep]).isdisjoint(used_prior_seeds)
        initial = [s for s in allstates if s['replicate'] == rep and s['mechanical_time_s'] == 0]
        assert len({s['canonical_geometry'] for s in initial}) == 1
    assert len(allstates) == NOMINAL_STATES and len(unique) == UNIQUE_STATES
    states = [s for s in allstates if s['unique_geometry']]
    jobs = []
    for st in states:
        sid, rep = st['state_id'], st['replicate']
        dest = OUT/'states'/sid
        poolpath = OUT/'shared'/f'release_candidates_rep{rep}.tsv'
        d = query(st['scene'], poolpath, dest/'candidate_clearance.tsv')
        assert len(d) == POOL
        ids = np.flatnonzero(d >= RADIUS+1e-12)[:N]
        assert len(ids) == N, 'Insufficient valid central starts: no fallback is permitted'
        xyz = np.loadtxt(poolpath)[ids]
        with np.load(dest/'geometry.npz', allow_pickle=False) as a:
            mask, extent, dx = a['common_support'], a['extent_um'], float(a['dx_um'])
        ij = np.floor((xyz[:, :2]-[extent[0], extent[2]])/dx).astype(int)
        assert mask[ij[:, 1], ij[:, 0]].all(), 'Central starts outside common support; no silent filtering'
        np.savetxt(st['initial'], xyz, fmt='%.17g')
        np.savez_compressed(dest/'initial_metadata.npz', xyz_um=xyz, candidate_id=ids, clearance_um=d[ids])
        for offset in range(0, N, 64):
            ip = OUT/'initials'/f'{sid}_{offset:03d}.tsv'
            np.savetxt(ip, xyz[offset:offset+64], fmt='%.17g')
            jobs.append(dict(name=f'{sid}_{offset:03d}', state_id=sid, method='reflect_dt15',
                scene=st['scene'], domain=st['domain'], initial=str(ip),
                output=str(OUT/'transport'/f'{sid}_{offset:03d}.csv'), N=64, T=HORIZON, dt=DT, seed=st['seed'], offset=offset))
        print(f'INITIALS: {sid}; 256 valid centers prepared.', flush=True)
    for st in allstates:
        if not st['unique_geometry']:
            source = OUT/'states'/st['canonical_geometry']
            shutil.copy2(source/'initial.tsv', st['initial'])
            shutil.copy2(source/'initial_metadata.npz', OUT/'states'/st['state_id']/'initial_metadata.npz')
    csvsave(OUT/'cohort_inventory.csv', inventory)
    save(OUT/'states.json', allstates)
    for p in list((OUT/'states').glob('*/*')) + list((OUT/'shared').glob('*')) + list((OUT/'initials').glob('*')):
        if p.is_file(): protect(p)
    protect(OUT/'cohort_inventory.csv')
    protect(OUT/'states.json')
    assert all(sha(p) == h for p, h in protected.items()), 'Input changed during preparation'
    protocol = dict(status='fixed allocation; pending transport, no outcome-dependent stopping',
        prior_protocol=str(OLD/'protocol.json'), replicas=list(REPLICAS), mechanical_layer_seeds=seedsets,
        states=states, jobs=jobs, nominal_states=NOMINAL_STATES, unique_states=UNIQUE_STATES, total_unique_histories=N*UNIQUE_STATES,
        workers=4, N=N, T=HORIZON, dt_s=DT, radius_um=RADIUS, D_um2_s=DIFFUSION, H_um=H, L_um=L, layers=50,
        mechanical_times_s=list(TIMES), f0_nN=[0, 2.2, 4.4, 8.8],
        release='same RNG rule [916002,replicate]; first256 valid centers from fixed65536 proposals in radius1.19um cylinder, z2..3um; not connectivity filtered',
        pairing='same proposal pool and per-index Brownian streams within replica; accepted coordinates can differ',
        geometry='continuous backbone/cohesion capsules and platelet disks extruded by H; native256x256 intersection of layer supports',
        boundaries='geometric reflection at solids and lateral boundary; absorbing z0 and z5; lateral mask not radius-eroded',
        connectivity=dict(subdivisions_xy=2, nz=100, hxy_um=47.6/512, hz_um=.05, method='inner6/outer26 plus continuous bridge witnesses; not converged connectivity'),
        geometry_MC=dict(N=MC, seed_rule='[919201,replica]', region='full47.6x47.6x5 box conditioned on common support'),
        failures='fatal, never censoring; partial outputs require inspection',
        temporal_equivalence='still inconclusive; finest tested step retained as descriptive pilot',
        old_replicas_not_rerun=True, mechanics_rerun=False, manuscript_changed=False, zip_changed=False,
        scientific_analysis_deferred=True, source_sha256=protected)
    save(OUT/'protocol.json', protocol)
    assert len(jobs) == TRANSPORT_JOBS
    status('prepared', unique_states=UNIQUE_STATES, transport_jobs=TRANSPORT_JOBS)
    print(f'PREPARED: replicas {REPLICAS}; {UNIQUE_STATES} unique geometries, {TRANSPORT_JOBS} chunks, {N*UNIQUE_STATES} histories; no mechanics rerun.', flush=True)


def geometry_job(s):
    dest = OUT/'states'/s['state_id']
    summary = conn.compute(OUT/'source_snapshot/common_connectivity.exe', Path(s['scene']), Path(s['domain']), Path(s['initial']), dest/'connectivity')
    points = OUT/'shared'/f"volume_rep{s['replicate']}.tsv"
    d = query(s['scene'], points, dest/'volume_clearance.tsv', .26)
    xyz = np.loadtxt(points)
    with np.load(dest/'geometry.npz', allow_pickle=False) as a:
        mask, extent, dx = a['common_support'], a['extent_um'], float(a['dx_um'])
    ij = np.floor((xyz[:, :2]-[extent[0], extent[2]])/dx).astype(int)
    inside = mask[ij[:, 1], ij[:, 0]]
    np.savez_compressed(dest/'volume_samples.npz', inside_common=inside, void=inside & (d > 0),
                        accessible=inside & (d >= RADIUS), proposal_ids=np.arange(MC), radius_um=RADIUS)
    return dict(state_id=s['state_id'], **summary)


def run():
    p = load(OUT/'protocol.json')
    assert all(sha(path) == h for path, h in p['source_sha256'].items())
    transport.EXE = OUT/'source_snapshot/transport_reflection_long.exe'
    started = time.time()
    done, geometry = [], []
    assert len(p['states']) == UNIQUE_STATES and len(p['jobs']) == TRANSPORT_JOBS
    status('running', phase='geometry_and_transport', geometry_done=0, geometry_total=UNIQUE_STATES, transport_done=0, transport_total=TRANSPORT_JOBS)
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(geometry_job, s): 'geometry' for s in p['states']}
        futures.update({pool.submit(transport.execute, j): 'transport' for j in p['jobs']})
        for f in as_completed(futures):
            kind = futures[f]
            try:
                result = f.result()
            except BaseException:
                for pending in futures: pending.cancel()
                raise
            (geometry if kind == 'geometry' else done).append(result)
            status('running', phase='geometry_and_transport', geometry_done=len(geometry), geometry_total=UNIQUE_STATES,
                   transport_done=len(done), transport_total=TRANSPORT_JOBS, elapsed_s=time.time()-started)
            print(f'COMPLETED {kind}: geometry {len(geometry)}/{UNIQUE_STATES}; transport {len(done)}/{TRANSPORT_JOBS}.', flush=True)
    assert all(sha(path) == h for path, h in p['source_sha256'].items())
    csvsave(OUT/'connectivity_summary.csv', geometry)
    save(OUT/'execution_summary.json', dict(completed=True, jobs=done, connectivity=geometry,
        wall_s=time.time()-started, source_unchanged=True, scientific_analysis_deferred=True))
    status('completed_awaiting_analysis', geometry_done=UNIQUE_STATES, transport_done=TRANSPORT_JOBS, wall_s=time.time()-started)
    print('FINISHED: raw transport and geometry available. No manuscript, report or ZIP updates.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['prepare', 'run', 'all'])
    args = parser.parse_args()
    try:
        if args.mode in ('prepare', 'all'): prepare()
        if args.mode in ('run', 'all'): run()
    except BaseException as exc:
        if OUT.exists(): status('failed', error=str(exc))
        raise
