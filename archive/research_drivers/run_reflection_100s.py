"""Fixed 100 s extension. Scientific analysis only after all jobs finish.

Same per-particle RNG streams as the 10 s pilot; four local CPU processes.
No mechanics rerun, no optional stopping, no changes to production Program19.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib, json, subprocess, sys, time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'outputs/verificacion_reflexion_20260919_v2'
OUT=ROOT/'outputs/verificacion_reflexion_100s_20260919'
CPP=ROOT/'work/transport_reflection_long.cpp'
EXE=OUT/'transport_reflection_long.exe'
METHODS={'reflect_dt31':.00003125,'reflect_dt15':.000015625}

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def load(p):return json.loads(p.read_text())
def data(p):return np.genfromtxt(p,delimiter=',',names=True,ndmin=1)
def execute(j):
    path=Path(j['output']);done=path.with_suffix('.done.json')
    if done.exists():
        q=load(done);assert q['job']==j and q['sha256']==sha(path)
        return q
    if path.exists():raise RuntimeError(f'Partial output requires inspection: {path}')
    args=[EXE,'run',j['scene'],j['domain'],j['initial'],path,j['N'],j['T'],j['dt'],j['seed'],.25,.4,0,j['offset']]
    start=time.perf_counter()
    with path.with_suffix('.log').open('w') as f:
        p=subprocess.run(list(map(str,args)),stdout=f,stderr=subprocess.STDOUT)
    q=dict(job=j,exit_code=p.returncode,wall_s=time.perf_counter()-start)
    if p.returncode:
        save(path.with_suffix('.failure.json'),q)
        raise RuntimeError(f"Execution failure; not censoring: {j['name']}")
    # Structural integrity only. Do not inspect scientific outcomes here.
    d=data(path)
    assert len(d)==j['N'] and np.array_equal(d['particle'],np.arange(j['offset'],j['offset']+j['N']))
    assert all(np.isfinite(d[n]).all() for n in d.dtype.names)
    assert np.isin(d['outcome'],[0,1,2]).all()
    assert ((d['time_s']>0)&(d['time_s']<=j['T'])).all()
    q['sha256']=sha(path);save(done,q)
    print(json.dumps(dict(job=j['name'],status='completed',wall_s=q['wall_s'])),flush=True)
    return q

def prepare():
    OUT.mkdir(exist_ok=True)
    for name in ['chunks','initials','regression']:(OUT/name).mkdir(exist_ok=True)
    controls=load(OLD/'controls.json');assert controls['passed']
    assert sha(OLD/'transport_reflection_v2.exe')==controls['exe_sha256']
    assert all(sha(Path(p))==h for p,h in controls['source_sha256'].items())
    source=(ROOT/'work/transport_reflection_v2.cpp').read_text()
    expected=source.replace('if(mode!="run"||argc!=13)', 'if(mode!="run"||(argc!=13&&argc!=14))')
    expected=expected.replace('seed radius D reject");','seed radius D reject [particle_offset]");')
    needle=' double r=std::stod(argv[10]),D=std::stod(argv[11]);bool rejection=std::stoi(argv[12])!=0;'
    expected=expected.replace(needle,needle+'\n int offset=argc==14?std::stoi(argv[13]):0;if(offset<0)throw std::runtime_error("negative particle offset");')
    expected=expected.replace('seed_hash(i))','seed_hash(i+offset))').replace('out<<i<<","<<time','out<<i+offset<<","<<time')
    assert CPP.read_text()==expected,'Only offset/RNG indexing changes are allowed'
    subprocess.run(['C:/Strawberry/c/bin/g++.exe','-std=c++17','-O3','-DNDEBUG','-static',str(CPP),'-o',str(EXE)],check=True)
    original=load(OLD/'verification_protocol.json')
    j0=original['jobs'][0];initial=np.loadtxt(j0['initial'])
    tests=[]
    for name,dt in METHODS.items():
        for offset in [0,32,128,224]:
            ip=OUT/'regression'/f'{name}_{offset}_initial.tsv';np.savetxt(ip,initial[offset:offset+4],fmt='%.17g')
            tests.append(dict(name=f'replay_{name}_{offset}',scene=j0['scene'],domain=j0['domain'],initial=str(ip),
                output=str(OUT/'regression'/f'{name}_{offset}.csv'),N=4,T=10,dt=dt,seed=919010,offset=offset))
    with ThreadPoolExecutor(max_workers=4) as pool:replays=list(pool.map(execute,tests))
    for q in replays:
        j=q['job'];name='reflect_dt31' if j['dt']==METHODS['reflect_dt31'] else 'reflect_dt15'
        old=data(OLD/'real_pilot'/f'{name}.csv')[j['offset']:j['offset']+4]
        new=data(j['output']);assert old.dtype==new.dtype
        for field in old.dtype.names:assert np.array_equal(old[field],new[field]),(j['name'],field)
    save(OUT/'regression.json',dict(passed=True,trajectories=32,comparison='All fields bit-identical to 10 s pilot, including global particle offsets',jobs=replays))
    jobs=[]
    # Fine jobs first to avoid a long fine-only tail with four workers.
    for name in ['reflect_dt15','reflect_dt31']:
        for offset in range(0,256,64):
            ip=OUT/'initials'/f'initial_{offset:03d}.tsv';np.savetxt(ip,initial[offset:offset+64],fmt='%.17g')
            jobs.append(dict(name=f'{name}_{offset:03d}',method=name,scene=j0['scene'],domain=j0['domain'],initial=str(ip),
                output=str(OUT/'chunks'/f'{name}_{offset:03d}.csv'),N=64,T=100,dt=METHODS[name],seed=919010,offset=offset))
    sources=[CPP,ROOT/'work/transport_geometry_bridge.cpp',EXE,Path(j0['scene']),Path(j0['domain']),Path(j0['initial']),Path(__file__)]
    sources+=sorted((OUT/'initials').glob('*.tsv'))
    protocol=dict(status='fixed before extended outcome inspection',N=256,horizon_s=100,max_workers=4,jobs=jobs,
        frozen_mechanical_time_s=10,layer_count=50,H_um=.1,L_um=5,radius_um=.25,D_um2_s=.4,
        primary_endpoints=['S(100 s) coarse-minus-fine','relative RMST(100 s) coarse-minus-fine'],
        descriptive_times_s=[10,25,50,75,100],margins=dict(S_absolute=.02,RMST_relative=.03),
        confidence='97.5% marginal CIs; Bonferroni nominal joint95% for two primary endpoints at100s only',
        S_interval='Conservative exact discordant-pair binomial bounds, alpha=.025 and tails alpha/4',
        RMST_interval='20000 paired bootstrap resamples; seed919031; quantiles .0125 and .9875',
        stopping_rule='Fixed N256 and100s; no outcome-dependent enlargement',
        inheritance='Same positions and per-particle fine Gaussian streams as10s pilot; extension is not independent replication',
        mechanics_rerun=False,production_changed=False,spatial_refinement=False,reject_comparator_rerun=False,
        geometry='Continuous extruded solid primitives incl cohesion; unchanged native raster lateral common domain',
        failure_policy='Explicit failure, never censored. Analysis blocked unless all jobs complete.',
        analysis_policy='No interim scientific analysis. Analyze once all eight chunks complete.',
        source_sha256={str(p):sha(p) for p in sources})
    p=OUT/'protocol.json'
    if p.exists():assert load(p)==protocol
    else:save(p,protocol)
    print('Prepared: fixed100s,256paired particles,two temporal steps,fourCPU processes;32replays exact.',flush=True)

def run():
    p=load(OUT/'protocol.json')
    assert all(sha(Path(k))==v for k,v in p['source_sha256'].items())
    start=time.time();save(OUT/'status.json',dict(status='running',started_unix=start,jobs_done=0,jobs_total=len(p['jobs'])))
    results=[]
    try:
        with ThreadPoolExecutor(max_workers=p['max_workers']) as pool:
            futures=[pool.submit(execute,j) for j in p['jobs']]
            for f in as_completed(futures):
                results.append(f.result())
                save(OUT/'status.json',dict(status='running',started_unix=start,jobs_done=len(results),jobs_total=len(p['jobs'])))
    except Exception as e:
        save(OUT/'status.json',dict(status='failed',error=str(e),jobs_done=len(results)));raise
    assert all(sha(Path(k))==v for k,v in p['source_sha256'].items())
    save(OUT/'execution_summary.json',dict(completed=True,wall_s=time.time()-start,jobs=results,sources_unchanged=True))
    save(OUT/'status.json',dict(status='all_simulations_completed',jobs_done=len(results),wall_s=time.time()-start))
    # Results are analyzed for the first time only here.
    import report_reflection_100s
    report_reflection_100s.main()
    save(OUT/'status.json',dict(status='completed_and_reported',jobs_done=len(results),wall_s=time.time()-start))

if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='run':run()
    else:raise ValueError('Use prepare or run')
