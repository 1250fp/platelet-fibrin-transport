"""Fixed-allocation verification jobs; at most two simultaneous processes."""
from pathlib import Path
import concurrent.futures,json,subprocess,time,hashlib
import numpy as np
from verify_reflection_v1 import ROOT,OUT,EXE,CPP,SOURCE,domain,save,sha

def job(j):
    path=Path(j['output']);path.parent.mkdir(exist_ok=True)
    done=path.with_suffix('.done.json')
    if done.exists():return json.loads(done.read_text())
    if path.exists():raise RuntimeError(f'Incomplete existing result: inspect {path}')
    args=[EXE,'run',j['scene'],j['domain'],j['initial'],path,j['N'],j['T'],j['dt'],j['seed'],j['radius'],j['D'],j['reject']]
    started=time.perf_counter()
    with path.with_suffix('.log').open('w') as log:
        result=subprocess.run(list(map(str,args)),stdout=log,stderr=subprocess.STDOUT)
    summary=dict(name=j['name'],exit_code=result.returncode,wall_s=time.perf_counter()-started,
                 output=str(path),job=j)
    if result.returncode:
        save(path.with_suffix('.failure.json'),summary)
        raise RuntimeError(f"{j['name']} failed: {path.with_suffix('.log')}")
    data=np.genfromtxt(path,delimiter=',',names=True)
    assert len(data)==j['N']
    summary.update(escaped=int((data['outcome']!=0).sum()),censored=int((data['outcome']==0).sum()),
                   sha256=sha(path))
    save(done,summary);print(json.dumps({k:v for k,v in summary.items() if k!='job'}),flush=True)
    return summary

def main():
    controls=json.loads((OUT/'controls.json').read_text());assert controls['passed']
    assert sha(EXE)==controls['exe_sha256']
    real=SOURCE/'states/rep1_f0_2.2_tm10'
    init=real/'initial_batch00.tsv';scene=real/'solid_full.tsv';reg=OUT/'real_domain.txt'
    (OUT/'real_pilot').mkdir(exist_ok=True);(OUT/'slab_control').mkdir(exist_ok=True)
    empty=OUT/'controls/empty.tsv'
    slabreg=OUT/'slab_control/domain.txt';slabinit=OUT/'slab_control/initial.tsv'
    domain(slabreg,np.ones((40,2),bool),-.05,-1,.05,1)
    np.savetxt(slabinit,np.tile([0,0,.5],(4096,1)),fmt='%.17g')
    jobs=[]
    for name,dt,reject in [('reflect_dt31',.00003125,0),('reflect_dt15',.000015625,0),('reject_dt15',.000015625,1)]:
        jobs.append(dict(name=name,scene=str(scene),domain=str(reg),initial=str(init),
            output=str(OUT/'real_pilot'/f'{name}.csv'),N=256,T=10,dt=dt,seed=919010,radius=.25,D=.4,reject=reject))
    for dt in [.0005,.000125,.00003125,.000015625]:
        name='slab_reflect_dt'+str(dt)
        jobs.append(dict(name=name,scene=str(empty),domain=str(slabreg),initial=str(slabinit),
            output=str(OUT/'slab_control'/f'{name}.csv'),N=4096,T=2,dt=dt,seed=919020,radius=.25,D=.4,reject=0))
    # Same first 4096 initial positions and noise streams for a reject comparator.
    jobs.append(dict(name='slab_reject_dt0.0005',scene=str(empty),domain=str(slabreg),initial=str(slabinit),
            output=str(OUT/'slab_control/slab_reject_dt0.0005.csv'),N=4096,T=2,dt=.0005,seed=919020,radius=.25,D=.4,reject=1))
    protocol=dict(status='fixed allocation before outcome inspection',jobs=jobs,max_workers=2,
        mechanics_rerun=False,production_changed=False,
        reflection='piecewise-linear increments with normal specular reflection; constant D, no drift; curved/multiple surfaces require temporal verification',
        real_primary=['S(10s) coarse-minus-fine difference','relative RMST(10s) coarse-minus-fine difference'],
        practical_margins=dict(S_absolute=.02,RMST_relative=.03),
        simultaneous_confidence='Bonferroni 95% over two primary endpoints: 97.5% marginal intervals',
        S_interval='conservative exact binomial bounds on both discordant-pair probabilities; not zero-width bootstrap',
        RMST_interval='paired percentile bootstrap ratio, 20000 resamples, 1.25% and 98.75% quantiles',
        stop_rule='256 paired histories at10s; no outcome-dependent sample enlargement; inconclusive is not passed',
        escape_clock='grid-end time is primary; linear segment crossing time saved separately, not asserted exact Brownian hitting time',
        numerical_failures='fatal errors, never reclassified as censoring',
        source_sha256={str(p):sha(p) for p in [scene,init,reg,EXE,CPP]})
    p=OUT/'verification_protocol.json'
    if p.exists():assert json.loads(p.read_text())==protocol
    else:save(p,protocol)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(job,j) for j in jobs]
        results=[f.result() for f in futures]
    assert all(sha(Path(p))==h for p,h in protocol['source_sha256'].items())
    save(OUT/'execution_summary.json',dict(completed=True,jobs=results,sources_unchanged=True))

if __name__=='__main__':main()
