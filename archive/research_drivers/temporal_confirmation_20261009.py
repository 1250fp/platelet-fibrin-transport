"""Fixed-allocation temporal confirmation; immutable input geometries, eight processes.

No mechanical simulation, interim outcome analysis, automatic enlargement, or paper edit.
Only the finest Gaussian interval differs from the production transport executable.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/temporal_confirmation_20261009'
OLD = ROOT / 'outputs/verificacion_reflexion_100s_20260919'
SNAP = OUT / 'source_snapshot'
EXE = SNAP / 'transport_temporal.exe'
PYTHON = sys.executable
N, WORKERS, CHUNK, HORIZON = 8192, 8, 64, 100.
STEPS = {'dt31': .00003125, 'dt15': .000015625, 'dt7': .0000078125}
SELECTION = [
    ('piloto_homogeneo_reflexion_20260919', 'rep1_f0_0_tm0', 'Initial state, before mechanical evolution'),
    ('piloto_homogeneo_reflexion_20260919', 'rep1_f0_2.2_tm10', 'Same geometry as the previous temporal pilot'),
    ('transporte_nueve_replicas_15s_20261001', 'rep9_f0_8.8_tm15', 'Latest mechanical state, largest baseline, another replica'),
]


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')
    tmp.replace(path)


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def status(stage, **extra):
    save(OUT/'status.json', dict(status=stage, pid=os.getpid(),
         updated_local=time.strftime('%Y-%m-%dT%H:%M:%S'), **extra))


def call(args, log):
    with Path(log).open('w', encoding='utf-8') as handle:
        subprocess.run(list(map(str, args)), stdout=handle, stderr=subprocess.STDOUT,
                       check=True, creationflags=subprocess.CREATE_NO_WINDOW)


def read_events(path):
    return np.genfromtxt(path, delimiter=',', names=True, ndmin=1)


def validate_events(path, job):
    data = read_events(path)
    assert len(data) == job['N']
    assert np.array_equal(data['particle'], np.arange(job['offset'], job['offset']+job['N']))
    assert all(np.isfinite(data[k]).all() for k in data.dtype.names)
    assert np.isin(data['outcome'], [0, 1, 2]).all()
    assert ((data['time_s'] > 0) & (data['time_s'] <= job['T'])).all()
    assert (data['time_s'][data['outcome'] == 0] == job['T']).all()
    assert (data['zf'][data['outcome'] == 1] == 0).all()
    assert (np.abs(data['zf'][data['outcome'] == 2] - job.get('L', 5.)) < 1e-9).all()
    assert np.array_equal(np.c_[data['x0'], data['y0'], data['z0']],
                          np.loadtxt(job['initial'], ndmin=2))
    assert (data['min_clearance_um'] >= .25-2e-10).all()
    return data


def args_for(exe, job):
    return [exe, 'run', job['scene'], job['domain'], job['initial'], job['output'],
            job['N'], job['T'], job['dt'], job['seed'], .25, .4, 0, job['offset']]


def copy_exact(source, target):
    source, target = Path(source), Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        assert sha(source) == sha(target), f'Existing snapshot differs: {target}'
    else:
        shutil.copy2(source, target)


def verify_hashes(mapping):
    for path, expected in mapping.items():
        assert sha(path) == expected, f'Input hash changed: {path}'


@contextmanager
def locked():
    import msvcrt
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT/'execution.lock').open('a+b') as handle:
        if handle.tell() == 0:
            handle.write(b'0')
            handle.flush()
        handle.seek(0)
        try:
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        except OSError as error:
            raise RuntimeError('Another process owns this experiment') from error
        try:
            yield
        finally:
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)


def controls():
    dest = OUT/'controls'
    dest.mkdir(exist_ok=True)
    legacy = SNAP/'transport_legacy_replay.exe'
    compiler = 'C:/Strawberry/c/bin/g++.exe'
    source = SNAP/'transport_reflection_temporal_20261009.cpp'
    call([compiler, '-std=c++17', '-O3', '-DNDEBUG', '-static', source, '-o', EXE], dest/'compile.log')
    call([compiler, '-std=c++17', '-O3', '-DNDEBUG', '-static',
          '-DTEMPORAL_NOISE_DT=0.000015625', source, '-o', legacy], dest/'compile_legacy.log')
    original = (SNAP/'transport_reflection_long.cpp').read_text()
    changed = source.read_text()
    oldline = ' const double noise_dt=0.000015625;int stride=int(std::llround(dt/noise_dt));'
    replacement = ''' // Only the finest Gaussian interval changes in this verification executable.
 // A legacy build with TEMPORAL_NOISE_DT=0.000015625 must replay old paths exactly.
#ifndef TEMPORAL_NOISE_DT
#define TEMPORAL_NOISE_DT 0.0000078125
#endif
 const double noise_dt=TEMPORAL_NOISE_DT;int stride=int(std::llround(dt/noise_dt));'''
    assert original.count(oldline) == 1 and changed == original.replace(oldline, replacement)
    prior = load(OLD/'protocol.json')
    verify_hashes(prior['source_sha256'])
    replays = []
    for method in ('reflect_dt15', 'reflect_dt31'):
        full = read_events(OLD/f'{method}.csv')
        exemplar = next(j for j in prior['jobs'] if j['method'] == method)
        for offset in (0, 128):
            initial = dest/f'replay_{method}_{offset}.tsv'
            np.savetxt(initial, np.c_[full['x0'], full['y0'], full['z0']][offset:offset+2], fmt='%.17g')
            job = dict(exemplar, initial=str(initial), output=str(dest/f'replay_{method}_{offset}.csv'),
                       N=2, offset=offset)
            call(args_for(legacy, job), dest/f'replay_{method}_{offset}.log')
            got = validate_events(job['output'], job)
            expected = full[offset:offset+2]
            assert all(np.array_equal(got[k], expected[k]) for k in got.dtype.names)
            replays.append(dict(method=method, offset=offset, count=2, all_fields_bit_identical=True))
    # Actual C++ Gaussian streams: same terminal Brownian increment at all levels,
    # correct free-space moments, and invariance under particle chunking.
    empty = dest/'empty.tsv'
    empty.write_text('', encoding='ascii')
    domain = dest/'free_domain.txt'
    domain.write_text('2 2 -10000 -10000 10000 20000\n1 1\n1 1\n', encoding='ascii')
    initial = dest/'free_initial.tsv'
    np.savetxt(initial, np.tile([0., 0., 10000.], (4096, 1)), fmt='%.17g')
    free = {}
    for name, dt in STEPS.items():
        job = dict(scene=str(empty), domain=str(domain), initial=str(initial),
                   output=str(dest/f'free_{name}.csv'), N=4096, T=.001, dt=dt,
                   seed=2026100901, offset=0, L=20000.)
        call(args_for(EXE, job), dest/f'free_{name}.log')
        d = validate_events(job['output'], job)
        free[name] = np.c_[d['xf']-d['x0'], d['yf']-d['y0'], d['zf']-d['z0']]
        expected_variance = 2*.4*.001
        assert (np.abs(free[name].mean(axis=0)) < 6*np.sqrt(expected_variance/4096)).all()
        assert (np.abs(free[name].var(axis=0, ddof=1)/expected_variance-1) < .15).all()
    maxerr = max(float(np.max(np.abs(x-free['dt7']))) for x in free.values())
    assert maxerr < 1e-9
    split = dict(job, N=8, offset=128, output=str(dest/'free_split.csv'))
    split['initial'] = str(dest/'free_split_initial.tsv')
    np.savetxt(split['initial'], np.tile([0., 0., 10000.], (8, 1)), fmt='%.17g')
    call(args_for(EXE, split), dest/'free_split.log')
    splitdata = validate_events(split['output'], split)
    fulldata = read_events(dest/'free_dt7.csv')[128:136]
    assert all(np.array_equal(splitdata[k], fulldata[k]) for k in splitdata.dtype.names)
    # Existing independent geometric controls remain unchanged, hash-verified.
    geometric = load(ROOT/'outputs/verificacion_reflexion_20260919_v2/controls.json')
    assert geometric['passed']
    result = dict(passed=True, engine_change='finest Gaussian interval only',
                  legacy_replays=replays, free_paths_per_level=4096,
                  free_nested_endpoint_max_error_um=maxerr, offset_replay_bit_identical=True,
                  inherited_geometric_cases=geometric['total_cases'])
    save(OUT/'controls.json', result)
    return result


def prepare():
    if (OUT/'protocol.json').exists():
        protocol = load(OUT/'protocol.json')
        verify_hashes(protocol['input_sha256'])
        return
    assert not (OUT/'transport').exists(), 'Inspect previous incomplete preparation first'
    status('preparing')
    for folder in ('source_snapshot', 'states'):
        (OUT/folder).mkdir(exist_ok=True)
    protected = {}
    for name in ('temporal_confirmation_20261009.py', 'report_temporal_confirmation_20261009.py',
                 'transport_reflection_temporal_20261009.cpp', 'transport_reflection_long.cpp',
                 'transport_geometry_bridge.cpp'):
        path = ROOT/'work'/name
        copy_exact(path, SNAP/name)
        protected[str(path)] = sha(path)
    clearance = ROOT/'outputs/transporte_nueve_replicas_15s_20261001/source_snapshot/common_clearance.exe'
    copy_exact(clearance, SNAP/'common_clearance.exe')
    control = controls()
    selected, jobs = [], []
    for cohort, sid, reason in SELECTION:
        path = ROOT/'outputs'/cohort/'protocol.json'
        prior = load(path)
        state = next(s for s in prior['states'] if s['state_id'] == sid)
        rep = state['replicate']
        dest = OUT/'states'/sid
        dest.mkdir(exist_ok=True)
        for key, filename in (('scene', 'solid_full.tsv'), ('domain', 'domain.txt')):
            source = Path(state[key])
            assert prior['source_sha256'][str(source)] == sha(source)
            protected[str(source)] = sha(source)
            copy_exact(source, dest/filename)
        rng = np.random.default_rng([20261009, rep, 11])
        count = 262144
        angle = rng.uniform(0, 2*np.pi, count)
        radial = 1.19*np.sqrt(rng.random(count))
        z = rng.uniform(2, 3, count)
        proposals = np.c_[radial*np.cos(angle), radial*np.sin(angle), z]
        np.savetxt(dest/'proposals.tsv', proposals, fmt='%.17g')
        call([SNAP/'common_clearance.exe', 'points', dest/'solid_full.tsv', dest/'proposals.tsv',
              dest/'clearance.tsv', .4], dest/'clearance.log')
        distance = np.loadtxt(dest/'clearance.tsv')
        header = (dest/'domain.txt').read_text().splitlines()[0].split()
        nx, ny = map(int, header[:2])
        xmin, ymin, dx, L = map(float, header[2:])
        mask = np.loadtxt(dest/'domain.txt', skiprows=1, dtype=int).reshape(ny, nx)
        ij = np.floor((proposals[:, :2]-[xmin, ymin])/dx).astype(int)
        assert ((ij >= 0) & (ij < [nx, ny])).all()
        allowed = mask[ij[:, 1], ij[:, 0]].astype(bool)
        ids = np.flatnonzero((distance >= .25+1e-12) & allowed)[:N]
        assert len(ids) == N, 'Not enough valid starts; no fallback or sample-size change'
        xyz = proposals[ids]
        np.savetxt(dest/'initial.tsv', xyz, fmt='%.17g')
        np.savez_compressed(dest/'initial_provenance.npz', candidate_id=ids, clearance_um=distance[ids], xyz_um=xyz)
        item = dict(state_id=sid, replicate=rep, f0_nN=state['f0_nN'],
                    mechanical_time_s=state['mechanical_time_s'], selection_reason=reason,
                    prior_protocol=str(path), original_scene=state['scene'], original_domain=state['domain'],
                    scene=str(dest/'solid_full.tsv'), domain=str(dest/'domain.txt'),
                    initial=str(dest/'initial.tsv'), seed=2026100900+rep, L=L)
        selected.append(item)
        for offset in range(0, N, CHUNK):
            initial = dest/f'initial_{offset:05d}.tsv'
            np.savetxt(initial, xyz[offset:offset+CHUNK], fmt='%.17g')
            for method in ('dt7', 'dt15', 'dt31'):
                name = f'{sid}_{method}_{offset:05d}'
                jobs.append(dict(name=name, state_id=sid, method=method, scene=item['scene'],
                    domain=item['domain'], initial=str(initial), output=str(OUT/'transport'/f'{name}.csv'),
                    N=CHUNK, T=HORIZON, dt=STEPS[method], seed=item['seed'], offset=offset, L=L))
        print(f'PREPARED {sid}: {N} valid central positions; no connectivity filtering', flush=True)
    for path in SNAP.rglob('*'):
        if path.is_file():
            protected[str(path)] = sha(path)
    for path in (OUT/'states').rglob('*'):
        if path.is_file():
            protected[str(path)] = sha(path)
    protocol = dict(created_local=time.strftime('%Y-%m-%dT%H:%M:%S'),
        status='fixed before new transport outcomes; independent confirmation sample',
        states=selected, jobs=jobs, N_per_geometry=N, steps_s=STEPS,
        geometries=len(selected), total_trajectories=len(selected)*N*3,
        horizon_s=100, workers=WORKERS, chunks=len(jobs), radius_um=.25, D_um2_s=.4,
        layers=50, H_um=.1, L_um=5.,
        primary_comparison='dt15 minus dt7 at 100 s in EACH selected geometry',
        secondary_comparison='dt31 minus dt15; diagnostic, not a requirement for accepting dt15',
        primary_endpoints=['S(100 s)', 'relative RMST(100 s) difference, finer-step denominator'],
        margins={'survival_absolute': .02, 'rmst_relative': .03},
        confidence='Same methods as historical reflection pilot: marginal 97.5% intervals; nominal joint95% within each endpoint pair. No simultaneous CI coverage claim across geometries.',
        survival_interval='Difference of exact Clopper-Pearson bounds for discordant-plus and discordant-minus probabilities; tail alpha/4, alpha=.025',
        rmst_interval='20000 paired bootstrap resamples; percentile .0125,.9875; fixed analysis seed2026100917',
        decision='dt15 acceptable on tested geometries only if both primary CIs lie strictly inside original margins in ALL three states; no extrapolation to all117 geometries or new stereological models',
        sample_size_rationale='Old 59/256 discordance predicts exact-CI halfwidth about5.08pp at1024,2.52pp at4096,1.77pp at8192 if zero bias.8192 fixed, not guaranteed power or equivalence.',
        descriptive_horizons_s=[10,25,50,75,100],
        initial_distribution='Independent RNG[20261009,replicate,11], uniform cylinder radius1.19um z2..3um conditioned on common lateral domain and solid distance>=0.25+1e-12. First8192 valid among262144 proposals.',
        brownian_coupling='Independent new transport seed2026100900+replicate; per-particle offset streams; coarse increments sum2 or4 finest increments. Same positions and Brownian path across all3steps.',
        unchanged='Existing continuous solids including cohesive links, native raster lateral support, reflecting solid/lateral boundaries, absorbing center at axial faces, end-of-step exit-time convention',
        stopping='Fixed allocation; no interim scientific analysis, no outcome-dependent enlargement or margin changes',
        failures='Explicit numerical or integrity failure, never censored; completed chunks retained; partial chunks require inspection',
        mechanics_rerun=False, old_transport_rerun=False, production_modified=False,
        manuscript_modified=False, zip_modified=False, connectivity_recalculated=False,
        controls=control, input_sha256=protected)
    save(OUT/'protocol.json', protocol)
    (OUT/'transport').mkdir()
    write_protocol_note(protocol)
    status('prepared', jobs_done=0, jobs_total=len(jobs))


def write_protocol_note(p):
    lines = ['# Prueba temporal confirmatoria del transporte', '',
        '## Pregunta y alcance', '',
        '¿El paso utilizado en el minipaper, 0.015625 ms, produce S(100 s) y RMST(100 s) prácticamente equivalentes a 0.0078125 ms en tres geometrías representativas ya guardadas?', '',
        'Es una verificación numérica del transporte congelado. No calibra el tiempo mecánico, la densidad de plaquetas, la representación 3D–2D ni la validez biológica de las capas.', '',
        '## Diseño fijado antes de ejecutar', '',
        '| Estado | Motivo |', '|---|---|']
    for s in p['states']:
        lines.append(f"| {s['state_id']} | {s['selection_reason']} |")
    lines += ['', '- Tres pasos: 0.03125, 0.015625 y 0.0078125 ms.',
        '- 8192 posiciones válidas por geometría; las mismas en los tres pasos. 73728 trayectorias en total.',
        '- Nuevas semillas para las posiciones y el ruido; no se reutilizan resultados exploratorios como si fueran confirmatorios.',
        '- Incrementos brownianos anidados: el paso doble suma dos incrementos del fino; el cuádruple suma cuatro.',
        '- Geometrías, radio 0.25 µm, D=0.4 µm²/s, 50 capas, H=100 nm y frontera lateral intactos.',
        '- Horizonte de transporte fijo: 100 s. Ocho procesos externos, sin paralelismo interno.',
        '- Ninguna evolución mecánica nueva ni modificación del número de plaquetas.', '',
        '## Criterios', '',
        r'$$\Delta S=S_{0.015625}(100)-S_{0.0078125}(100).$$', '',
        r'$$\delta_R=\frac{\mathrm{RMST}_{0.015625}(100)-\mathrm{RMST}_{0.0078125}(100)}{\mathrm{RMST}_{0.0078125}(100)}.$$', '',
        'Los intervalos completos deben quedar dentro de ±2 puntos porcentuales y ±3%, respectivamente. Se conservan los IC marginales del97.5% del piloto previo: cotas binomiales exactas conservadoras para la discordancia y20000 remuestreos pareados para RMST. Cobertura conjunta nominal95% para los dos indicadores de cada comparación; no se afirma cobertura simultánea de todos los IC entre geometrías.', '',
        'El paso actual solo se considerará apoyado en este conjunto si ambos criterios primarios pasan en cada una de las tres geometrías. El contraste0.03125–0.015625 es diagnóstico. No se exige que un paso más grueso pase para aceptar el paso intermedio. Las curvas S(t) serán descriptivas, no una prueba de equivalencia uniforme.', '',
        '## Por qué no solo1024', '',
        'La discordancia observada en el piloto terminado fue59/256. Manteniendo el método conservador de IC y suponiendo diferencia real nula, la semianchura orientativa es5.08pp para1024,2.52pp para4096 y1.77pp para8192. Por eso se fija8192; esto no asegura aprobación si existe sesgo o mayor discordancia.', '',
        '## Decisión y límites', '',
        '- **Cumple:** respaldar0.015625ms solo para las geometrías, radio, difusión, fronteras y horizonte comprobados.',
        '- **Inconcluso:** distinguir falta de precisión de una diferencia relevante; no aumentar automáticamente muestra ni márgenes.',
        '- **Fuera del margen:** revisar el integrador/refinar el paso; no actualizar resultados de producción automáticamente.',
        '- Un fallo numérico detiene el análisis y nunca cuenta como censura.',
        '- Las incertidumbres de estereología y conectividad permanecen pendientes por separado.', '',
        '## Trazabilidad', '',
        '`protocol.json` fija posiciones, semillas, parámetros, criterios y hashes. `controls.json` recoge replays bit a bit del ejecutor anterior y controles del acoplamiento browniano. `status.json` informa solo del estado operativo. Los resultados se analizan una sola vez al completar todos los bloques.', '',
        'Fuentes metodológicas: [Lakens2017](https://doi.org/10.1177/1948550617697177), sobre equivalencia y márgenes prefijados; [Giles, métodos Monte Carlo multinivel](https://people.maths.ox.ac.uk/gilesm/mlmc.html), sobre acoplar resoluciones mediante incrementos brownianos compartidos. Aquí se utiliza ese acoplamiento para comparar pasos; no se implementa un estimador multinivel.', '']
    (OUT/'protocolo.md').write_text('\n'.join(lines), encoding='utf-8')


def execute(job):
    output = Path(job['output'])
    marker = output.with_suffix('.done.json')
    if marker.exists():
        value = load(marker)
        assert value['job'] == job and value['sha256'] == sha(output)
        validate_events(output, job)
        return value
    if output.exists():
        raise RuntimeError(f'Partial output requires inspection; will not overwrite: {output}')
    started = time.perf_counter()
    try:
        call(args_for(EXE, job), output.with_suffix('.log'))
        validate_events(output, job)
        value = dict(job=job, exit_code=0, wall_s=time.perf_counter()-started, sha256=sha(output))
        save(marker, value)
        return value
    except Exception as error:
        save(output.with_suffix('.failure.json'), dict(job=job, error=str(error),
             wall_s=time.perf_counter()-started, numerical_failure_not_censoring=True))
        raise


def run():
    p = load(OUT/'protocol.json')
    verify_hashes(p['input_sha256'])
    complete, pending = [], []
    for job in p['jobs']:
        marker = Path(job['output']).with_suffix('.done.json')
        if marker.exists():
            complete.append(execute(job))
        else:
            assert not Path(job['output']).exists(), f"Partial chunk requires inspection: {job['output']}"
            pending.append(job)
    started = time.time()
    status('running', jobs_done=len(complete), jobs_total=len(p['jobs']), started_unix=started)
    try:
        iterator = iter(pending)
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            active = {}
            for _ in range(WORKERS):
                j = next(iterator, None)
                if j is not None:
                    active[pool.submit(execute, j)] = j
            while active:
                ready, _ = wait(active, return_when=FIRST_COMPLETED)
                failed = None
                for future in ready:
                    j = active.pop(future)
                    try:
                        result = future.result()
                        complete.append(result)
                        print(f"COMPLETED {len(complete)}/{len(p['jobs'])}: {j['name']} ({result['wall_s']:.1f}s)", flush=True)
                    except Exception as error:
                        failed = error
                if failed is not None:
                    status('failed', error=str(failed), jobs_done=len(complete), jobs_total=len(p['jobs']))
                    raise failed  # Already-running children finish; no new chunks submitted.
                status('running', jobs_done=len(complete), jobs_total=len(p['jobs']), started_unix=started,
                       active_jobs=[j['name'] for j in active.values()])
                for _ in range(WORKERS-len(active)):
                    j = next(iterator, None)
                    if j is not None:
                        active[pool.submit(execute, j)] = j
        verify_hashes(p['input_sha256'])
        save(OUT/'execution_summary.json', dict(completed=True, wall_s_this_session=time.time()-started,
             jobs=complete, jobs_total=len(p['jobs']), source_hashes_unchanged=True))
        status('all_simulations_completed', jobs_done=len(complete), jobs_total=len(p['jobs']))
        import report_temporal_confirmation_20261009 as report
        report.main()
        status('completed_and_reported', jobs_done=len(complete), jobs_total=len(p['jobs']))
    except Exception as error:
        status('failed', error=str(error), jobs_done=len(complete), jobs_total=len(p['jobs']))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['prepare', 'run', 'all'])
    options = parser.parse_args()
    with locked():
        if options.mode in ('prepare', 'all'):
            prepare()
        if options.mode in ('run', 'all'):
            run()
