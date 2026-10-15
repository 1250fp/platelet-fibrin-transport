"""Continue the nine paired mechanical replicas; no transport or manuscript edits."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['MPLBACKEND'] = 'Agg'
sys.dont_write_bytecode = True

PROGRAM = Path('external_data/program19_3')
EXPERIMENTS = PROGRAM / 'results/exploratory_experiments'
WORKSPACE = Path('external_data/research_workspace')
PREFIX = 'mechanics_10to15_R1_R9_'
FORCES = [0., 2.2, 4.4, 8.8]
WORKERS = 8


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding='utf-8')
    temp.replace(path)


def now():
    return datetime.now().isoformat(timespec='seconds')


def layers_of(run):
    manifest = run / 'stack_manifest.json'
    if manifest.exists():
        return [Path(p).resolve() for p in read(manifest)['layer_directories']]
    return sorted((run / 'mechanical_layers').glob('layer_*'))


def source_runs():
    sources = {}
    extension = EXPERIMENTS / 'f0_sweep_20260903_112342/paired_extension_to_50_layers/runs'
    for manifest in extension.glob('*/*/program19_3_extended_stack_*/stack_manifest.json'):
        run = manifest.parent
        config = read(run / 'config.json')
        rep = {4000: 1, 4020: 2}[config['seed']]
        key = (rep, float(config['filopodia_f0_nN']))
        assert key not in sources, key
        assert read(run / 'status.json')['status'] == 'completed'
        sources[key] = run
    for rep in range(3, 10):
        batches = list(EXPERIMENTS.glob(f'replica{rep:02d}_mechanics_*'))
        assert len(batches) == 1, (rep, batches)
        batch = batches[0]
        assert read(batch / 'progress.json')['status'] == 'completed'
        for row in read(batch / 'completed_runs.json'):
            assert row['status'] == 'completed'
            key = (rep, float(row['f0_nN']))
            assert key not in sources
            sources[key] = Path(row['run_directory'])
    assert set(sources) == {(rep, force) for rep in range(1, 10) for force in FORCES}
    return sources


def prepare():
    import numpy as np
    if list(EXPERIMENTS.glob(PREFIX + '*')):
        raise RuntimeError('A continuation batch already exists; inspect it instead of duplicating it.')
    sys.path.insert(0, str(PROGRAM))
    from program19_core.config import Program19Config
    from program19_core.persistence import load_state, load_rng
    from program19_core.state import state_arrays

    reference = EXPERIMENTS / ('f0_sweep_20260903_112342/runs/'
        'program19_3_20mechanical_layers_paper_consistent_thesis_audited_corrected_adaptive_seed4000_20260903_115152')
    old_hashes = {k.replace('\\', '/'): v for k, v in read(reference / 'manifest.json')['source_sha256'].items()}
    for name in ('config', 'profiles', 'geometry', 'state', 'wlc', 'bending', 'interactions',
                 'cohesion', 'integrator', 'mechanical_runner', 'parallel', 'persistence', 'layers'):
        rel = f'program19_core/{name}.py'
        assert sha(PROGRAM / rel) == old_hashes[rel], f'Mechanical code changed: {rel}'

    old_states = read(WORKSPACE / 'outputs/analisis_homogeneo_con_cohesion_20260917/states.json')
    old_layers = {(int(s['replicate']), float(s['f0_nN'])):
        [Path(p['path']).resolve() for p in s['layers_source']]
        for s in old_states if s['mechanical_time_s'] == 10}
    jobs = []
    protected = {}
    seed_blocks = {}
    nonphysical = {'seed', 'mechanical_layers', 'mechanical_workers', 'layer_spacing_um',
                   'output_root', 'run_transport'}
    for (rep, f0), source in sorted(source_runs().items()):
        cfg = read(source / 'config.json')
        layers = layers_of(source)
        assert len(layers) == cfg['mechanical_layers'] == 50
        assert cfg['duration_s'] == 10 and cfg['output_interval_s'] == .1
        assert cfg['save_trajectory'] and cfg['trajectory_detail'] == 'topology'
        assert cfg['save_intermediate_masks']
        if rep <= 2:
            assert layers == old_layers[(rep, f0)], 'Original cohort layer ordering changed'
        seeds = []
        for layer in layers:
            lc = read(layer / 'config.json')
            differences = {key for key in cfg if cfg[key] != lc[key]} - nonphysical
            assert not differences, (layer, differences)
            summary = read(layer / 'layer_summary.json')
            assert summary['status'] == 'completed'
            assert summary['seed'] == lc['seed']
            seeds.append(int(lc['seed']))
            state = load_state(layer / 'final_state.npz')
            assert math.isclose(state.time_s, 10., abs_tol=1e-10)
            with np.load(layer / 'checkpoint.npz', allow_pickle=False) as checkpoint:
                for key, value in state_arrays(state).items():
                    assert np.array_equal(value, checkpoint[key]), (layer, key)
            _, policy = load_rng(layer / 'checkpoint.npz', fallback_seed=np.random.SeedSequence(0))
            assert policy == 'exact_saved_rng_state'
            assert (layer / 'trajectory.npz').is_file()
            for name in ('config.json', 'final_state.npz', 'checkpoint.npz', 'layer_summary.json'):
                protected[str(layer / name)] = sha(layer / name)
        protected[str(source / 'config.json')] = sha(source / 'config.json')
        if (source / 'stack_manifest.json').exists():
            protected[str(source / 'stack_manifest.json')] = sha(source / 'stack_manifest.json')
        assert len(set(seeds)) == 50
        seed_blocks.setdefault(rep, seeds)
        assert seeds == seed_blocks[rep], 'Seed pairing changed'
        slug = f'rep{rep:02d}_f0_{f0:g}'.replace('.', 'p')
        jobs.append(dict(job_id=slug, replica=rep, f0_nN=f0, original_run=str(source),
            source_layers=list(map(str, layers)), seeds=seeds, source_config=cfg))
        print(f'{now()} VERIFIED {slug}: 50 states at 10 s and exact RNG checkpoints', flush=True)
    assert len({s for seeds in seed_blocks.values() for s in seeds}) == 450

    batch = EXPERIMENTS / (PREFIX + datetime.now().strftime('%Y%m%d_%H%M%S'))
    batch.mkdir(exist_ok=False)
    snapshot = batch / 'source_snapshot'
    snapshot.mkdir()
    shutil.copytree(PROGRAM / 'program19_core', snapshot / 'program19_core',
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copy2(__file__, snapshot / Path(__file__).name)
    for job in jobs:
        source_view = batch / 'source_views' / job['job_id']
        original = job.pop('source_config')
        adjusted = dict(original, layer_spacing_um=.1, mechanical_workers=WORKERS, run_transport=False)
        config = Program19Config.from_dict(adjusted)
        config.validate()
        save(source_view / 'config.json', config.to_dict())
        save(source_view / 'stack_manifest.json', dict(layer_directories=job['source_layers'],
            parent_run_directory=job['original_run'], operation='read-only continuation source index'))
        job['source_view'] = str(source_view)
        job['configuration_overrides'] = {k: {'old': original[k], 'new': adjusted[k]}
            for k in adjusted if original[k] != adjusted[k]}
        job['output_root'] = str(batch / 'runs' / job['job_id'])
    plan = dict(created_local=now(), operation='continue all nine paired replicas from 10 to 15 seconds',
        workers=WORKERS, stacks=36, layers_per_stack=50, total_layer_evolutions=1800,
        source_time_s=10, additional_time_s=5, target_time_s=15, output_interval_s=.1,
        H_um=.1, Lz_um=5., run_transport=False, jobs=jobs,
        source_sha256=protected,
        snapshot_sha256={str(p.relative_to(snapshot)).replace('\\', '/'): sha(p)
                         for p in snapshot.rglob('*.py')},
        source_view_sha256={str(p): sha(p) for p in (batch / 'source_views').rglob('*.json')},
        tests_passed=['checkpoint state/RNG round trip', 'serial/parallel continuation equivalence',
                      'completed-run continuation and catalogue'],
        notes=['Original results stay unchanged; child trajectories cover 10 to 15 s only.',
               'Original R1/R2 metadata used H=70nm; exported H=100nm matches the adopted analysis convention and does not change independent 2D mechanics.',
               'Automatic raster pore summaries are legacy diagnostics, not the continuous cohesion-aware paper transport analysis.'])
    save(batch / 'batch_manifest.json', plan)
    save(batch / 'progress.json', dict(status='prepared', completed_stacks=0, total_stacks=36,
        workers=WORKERS, updated_local=now(), run_transport=False))
    print('PREPARED_BATCH=' + str(batch), flush=True)


def run(batch):
    import numpy as np
    batch = Path(batch).resolve()
    with (batch / 'started.lock').open('x', encoding='utf-8') as handle:
        handle.write(f'pid={os.getpid()}\nstarted={now()}\n')
    plan = read(batch / 'batch_manifest.json')
    for rel, digest in plan['snapshot_sha256'].items():
        assert sha(batch / 'source_snapshot' / rel) == digest, rel
    for path, digest in plan['source_view_sha256'].items():
        assert sha(path) == digest, path
    sys.path.insert(0, str(batch / 'source_snapshot'))
    from program19_core.continuation import continue_mechanical_run
    current = dict(status='running', pid=os.getpid(), started_local=now(), updated_local=now(),
        completed_stacks=0, total_stacks=36, completed_layers_total=0, total_layers=1800,
        workers=WORKERS, run_transport=False)
    save(batch / 'progress.json', current)
    completed = []
    start = time.monotonic()
    try:
        for index, job in enumerate(plan['jobs']):
            relevant = [p for p in plan['source_sha256'] if p == str(Path(job['original_run']) / 'config.json')
                or p == str(Path(job['original_run']) / 'stack_manifest.json')
                or str(Path(p).parent) in job['source_layers']]
            for path in relevant:
                assert sha(path) == plan['source_sha256'][path], 'Source changed: ' + path
            root = Path(job['output_root'])
            assert not root.exists(), 'Existing output must be inspected, never silently repeated'
            current.update(job_id=job['job_id'], replica=job['replica'], f0_nN=job['f0_nN'],
                completed_layers_in_stack=0, active_workers=0, stage='starting',
                stack_progress_percent=0., completed_layers_total=index * 50, updated_local=now())
            save(batch / 'progress.json', current)
            print(f'{now()} START {index+1}/36 {job["job_id"]}: 10 -> 15 s; workers=8; transport=OFF', flush=True)
            last = [0., -1]

            def progress(payload):
                ts = time.monotonic()
                count = int(payload.get('completed_layers', current['completed_layers_in_stack']))
                if ts-last[0] < 60 and count == last[1] and payload.get('stage') != 'complete':
                    return
                current.update(stage=payload.get('stage'), updated_local=now(),
                    elapsed_s=round(ts-start, 1), completed_layers_in_stack=count,
                    completed_layers_total=index*50+count,
                    stack_progress_percent=round(100*float(payload.get('fraction', 0)), 2))
                for key in ('active_workers', 'parallel_workers', 'simulated_time_s'):
                    if key in payload:
                        current[key] = payload[key]
                save(batch / 'progress.json', current)
                print(f'{now()} {index+1}/36 {job["job_id"]}: {current["stage"]}; saved {count}/50', flush=True)
                last[:] = [ts, count]

            result = continue_mechanical_run(job['source_view'], additional_time_s=5.,
                output_root=root, mechanical_workers=WORKERS, progress_callback=progress)
            assert read(result / 'status.json')['status'] == 'completed'
            manifest = read(result / 'continuation_manifest.json')
            assert len(manifest['rng_policy_by_layer']) == 50
            assert all(x['policy'] == 'exact_saved_rng_state' for x in manifest['rng_policy_by_layer'])
            children = sorted((result / 'mechanical_layers').glob('layer_*'))
            assert len(children) == 50
            for parent, child, seed in zip(map(Path, job['source_layers']), children, job['seeds']):
                summary = read(child / 'layer_summary.json')
                assert summary['status'] == 'completed' and summary['seed'] == seed
                assert math.isclose(summary['simulated_time_s'], 15, abs_tol=1e-9)
                with np.load(parent / 'final_state.npz') as old, np.load(child / 'initial_state.npz') as initial:
                    for key in old.files:
                        assert np.array_equal(old[key], initial[key]), (child, key)
                with np.load(child / 'final_state.npz') as final, np.load(parent / 'final_state.npz') as old:
                    n = len(old['cohesion_pairs'])
                    for key in ('cohesion_pairs', 'cohesion_reference_lengths_um', 'cohesion_birth_times_s'):
                        assert np.array_equal(final[key][:n], old[key]), (child, key)
                assert not (child / 'first_passage.npz').exists()
            for path in relevant:
                assert sha(path) == plan['source_sha256'][path], 'Source changed during continuation: ' + path
            save(result / 'cohort_continuation_provenance.json', dict(original_run=job['original_run'],
                source_view=job['source_view'], replica=job['replica'], f0_nN=job['f0_nN'],
                source_layers=job['source_layers'], layer_seeds=job['seeds'],
                time_segments_s=[[0, 10], [10, 15]], H_um=.1, run_transport=False))
            completed.append(dict(replica=job['replica'], f0_nN=job['f0_nN'],
                source_run=job['original_run'], run_directory=str(result), status='completed',
                transport_run=False, source_time_s=10, final_time_s=15))
            save(batch / 'completed_runs.json', completed)
            current.update(completed_stacks=index+1, completed_layers_total=(index+1)*50)
            print(f'{now()} DONE {index+1}/36: {result}', flush=True)
        current.update(status='completed', stage='all_mechanics_complete', active_workers=0,
            completed_layers_in_stack=50, updated_local=now(), elapsed_s=round(time.monotonic()-start, 1))
        save(batch / 'execution_summary.json', dict(completed=True, stacks=36, layers=1800,
            target_time_s=15, workers=8, run_transport=False, elapsed_s=current['elapsed_s']))
        print(f'{now()} FINISHED: 36 stacks at 15 s; transport NOT run.', flush=True)
    except BaseException as exc:
        current.update(status='failed', error=str(exc), traceback=traceback.format_exc(), updated_local=now())
        raise
    finally:
        save(batch / 'progress.json', current)


if __name__ == '__main__':
    import multiprocessing
    multiprocessing.freeze_support()
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prepare', action='store_true')
    group.add_argument('--run', type=Path)
    args = parser.parse_args()
    if args.prepare:
        prepare()
    else:
        run(args.run)
