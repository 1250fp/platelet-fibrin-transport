"""Prepare matched frozen states; geometry only, never run transport implicitly.

All source runs and historical verification results are read-only. Outputs are
isolated. The local candidate preserves the verified ROI; global contraction is
reported separately. No cell graph is substituted for continuous transport.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.ndimage import generate_binary_structure, label

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / 'outputs/auditoria_trazabilidad_20260916'
OUT = ROOT / 'outputs/analisis_homogeneo_20260916'
CPP = ROOT / 'work/homogeneous_clearance.cpp'
EXE = CPP.with_suffix('.exe')
OLD_CPP = ROOT / 'work/transport_geometry_bridge.cpp'
TRANSPORT_EXE = ROOT / 'work/transport_continuous_refinement_dt15.exe'
H, HALF, L, RADIUS, D = .1, 2.975, 5., .25, .4
EPS, CAP = 1e-12, .4
SIDE, NZ = 128, 200  # Bounded bracketing graph, NOT a converged transport mask.
MC_N, CANDIDATE_N, N_INITIAL, BATCH_N = 65536, 65536, 1024, 256
protected = {}


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def protect(path):
    path = Path(path)
    protected[str(path)] = digest(path)
    return path


def read(path):
    return json.loads(protect(path).read_text(encoding='utf-8-sig'))


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')


def csvsave(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)


def shell(args):
    subprocess.run(list(map(str, args)), check=True, capture_output=True, text=True)


def points(scene, path, output, check=False):
    shell([EXE, 'check' if check else 'points', scene, path, output, CAP])
    return np.loadtxt(output)


def wilson(k, n):
    z = 1.959963984540054
    p = k/n; den = 1+z*z/n
    c = (p+z*z/(2*n))/den
    h = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0., c-h), min(1., c+h)]


def exit_labels(mask, connectivity):
    labels, n = label(mask, structure=generate_binary_structure(3, connectivity))
    exits = np.zeros(n+1, bool)
    exits[np.unique(labels[0])] = True; exits[np.unique(labels[-1])] = True
    exits[0] = False
    return labels, exits


def classification(dist, xyz):
    spacing = np.array([L/NZ, 2*HALF/SIDE, 2*HALF/SIDE])
    halfdiag = .5*np.linalg.norm(spacing)
    # Signed rounding guards bias both bounds outward, never toward a claim.
    safe = dist-halfdiag >= RADIUS+EPS
    possible = dist+halfdiag >= RADIUS-EPS
    assert np.all(~safe | possible)
    ij = np.floor((xyz[:, [2,1,0]]+[0,HALF,HALF])/spacing).astype(int)
    low_labels, low_exit = exit_labels(safe, 1)
    lower = low_exit[low_labels[tuple(ij.T)]]
    del low_labels
    upper_labels, upper_exit = exit_labels(possible, 3)
    upper = upper_exit[upper_labels[tuple(ij.T)]]
    assert np.all(~lower | upper)
    assert np.all(possible[tuple(ij.T)]), 'Continuously valid start excluded by outer mask'
    # 1: certified exit-connected; -1: disconnected within this ROI; 0: unresolved.
    status = np.where(lower, 1, np.where(~upper, -1, 0)).astype(np.int8)
    return status, dict(grid_side=SIDE, grid_z_cells=NZ, hxy_um=spacing[1], hz_um=spacing[0],
        half_cell_diagonal_um=halfdiag, safe_cell_volume_fraction=float(safe.mean()),
        possible_cell_volume_fraction=float(possible.mean()),
        F_conn_lower=float(lower.mean()), F_conn_upper=float(upper.mean()),
        connected_count=int(lower.sum()), disconnected_within_ROI_count=int((~upper).sum()),
        unresolved_count=int((upper & ~lower).sum()))


def tests():
    testdir = OUT/'tests'; testdir.mkdir(parents=True, exist_ok=True)
    shell(['C:/Strawberry/c/bin/g++.exe', '-std=c++17', '-O3', '-DNDEBUG', '-static', '-o', EXE, CPP])
    scene = testdir/'disk.tsv'
    np.savetxt(scene, [[0,0,0,0,.2,.5,1.5]], fmt='%.17g')
    rng = np.random.default_rng(916000)
    xyz = rng.uniform([-1,-1,0], [1,1,2], (2048,3))
    xyz = np.vstack([xyz, [0,0,1], [.3,0,1], [0,0,.4], [0,0,.5]])
    np.savetxt(testdir/'xyz.tsv', xyz, fmt='%.17g')
    got = points(scene, testdir/'xyz.tsv', testdir/'distance.tsv', True)
    analytic = np.minimum(CAP, np.hypot(np.maximum(np.linalg.norm(xyz[:,:2],axis=1)-.2,0),
        np.maximum(np.maximum(.5-xyz[:,2],xyz[:,2]-1.5),0)))
    assert np.allclose(got[:,0], analytic, rtol=0, atol=2e-15)
    assert np.allclose(got[:,0], got[:,1], rtol=0, atol=2e-15)
    # Independent vectorized capsule expression and original unaccelerated query.
    objects = rng.uniform(-1,1,(31,7)); objects[:,4]=rng.uniform(.02,.3,31)
    objects[:,5]=rng.uniform(0,1,31);objects[:,6]=objects[:,5]+rng.uniform(.01,.7,31)
    np.savetxt(testdir/'capsules.tsv',objects,fmt='%.17g')
    got = points(testdir/'capsules.tsv',testdir/'xyz.tsv',testdir/'capsule_distance.tsv',True)
    truth=np.full(len(xyz),CAP)
    for ax,ay,bx,by,r,zlo,zhi in objects:
        v=np.array([bx-ax,by-ay]);w=xyz[:,:2]-[ax,ay]
        u=np.clip((w@v)/np.dot(v,v),0,1)
        dxy=np.maximum(np.linalg.norm(w-u[:,None]*v,axis=1)-r,0)
        dz=np.maximum(np.maximum(zlo-xyz[:,2],xyz[:,2]-zhi),0)
        truth=np.minimum(truth,np.hypot(dxy,dz))
    assert np.allclose(got[:,0],truth,rtol=0,atol=2e-14)
    assert np.allclose(got[:,0],got[:,1],rtol=0,atol=2e-14)
    # Six-neighbor inner and 26-neighbor outer connectivity on exact controls.
    a=np.ones((8,9,9),bool);la,ex=exit_labels(a,1);assert ex[la].all()
    a[2]=False;a[5]=False;la,ex=exit_labels(a,3);assert not ex[la[3:5]].any()
    a[2,4,4]=True;la,ex=exit_labels(a,1);assert ex[la[3,4,4]]
    save(testdir/'results.json',dict(passed=True,analytic_disk_queries=len(xyz),
        random_capsule_queries=len(xyz),connectivity_controls=3,
        source_sha256=digest(CPP),included_source_sha256=digest(OLD_CPP),exe_sha256=digest(EXE)))


def prepare():
    assert not (OUT/'preparation_manifest.json').exists(), 'Preparation complete: do not overwrite'
    assert not OUT.exists(), 'Partial directory exists; inspect before retry'
    OUT.mkdir(parents=True)
    tick=time.perf_counter(); inv=read(AUDIT/'inventario.json'); layers=read(AUDIT/'capas_y_estados.json')
    assert not inv['issues']; grouped=defaultdict(list)
    for row in layers: grouped[row['stack_id']].append(row)
    protect(OLD_CPP); protect(TRANSPORT_EXE)
    test=read(ROOT/'outputs/transport_continuous_dt15_H100_20260911/tests.json')
    assert test['exe_sha256']==digest(TRANSPORT_EXE)
    # Freeze the exact executor and all input source identities in the proposal.
    for name in ['transport_continuous_refinement_dt15.cpp', 'transport_continuous_refinement_dt15.py']:
        protect(ROOT/'work'/name)
    globalrows=[]; states=[]
    for stack in inv['stacks']:
        sid=stack['id']; selected=sorted(grouped[sid],key=lambda q:q['layer'])
        assert len(selected)==50
        buffers={t:[] for t in (0.,5.,10.)}; metadata={t:[] for t in buffers}; regions={t:[] for t in buffers}
        for row in selected:
            layer=Path(row['path']); cfg=read(layer/'config.json')
            with np.load(protect(layer/'trajectory.npz'),allow_pickle=False) as f:
                nd,pl,pairs,offsets=f['nodes_um'],f['platelets_um'],f['backbone_pairs'],f['cohesion_frame_offsets']
                c_pairs=f['cohesion_pairs_flat']
                for snap in row['states']:
                    t=snap['target_time_s']; j=snap['frame_index']; k=row['layer']-1
                    a=nd[j,pairs[:,0]]; b=nd[j,pairs[:,1]]; p=pl[j]
                    fiber=np.c_[a,b,np.full(len(a),cfg['fiber_radius_um']),np.full(len(a),k*H),np.full(len(a),(k+1)*H)]
                    platelets=np.c_[p,p,np.full(len(p),cfg['platelet_radius_um']),np.full(len(p),k*H),np.full(len(p),(k+1)*H)]
                    buffers[t].append(np.vstack([fiber,platelets]))
                    metadata[t].append(dict(layer=row['layer'],seed=row['seed'],path=str(layer),
                        frame_index=j,saved_time_s=snap['saved_time_s'],backbone_primitives=len(a),platelets=len(p),
                        cohesive_links_saved=int(offsets[j+1]-offsets[j]),cohesion_is_geometric_obstacle=False))
            with np.load(protect(layer/'intermediate_masks.npz'),allow_pickle=False) as f:
                masks=f['analyzed_region_masks']
                for snap in row['states']: regions[snap['target_time_s']].append(masks[snap['frame_index']])
        for st in stack['states']:
            t=st['mechanical_time_s']; state_id=f'{sid}_tm{t:g}'; dest=OUT/'states'/state_id;dest.mkdir(parents=True)
            scene=np.concatenate(buffers[t]); assert np.isfinite(scene).all()
            np.savetxt(dest/'solid_full.tsv',scene,fmt='%.17g')
            # Keep the halo wider than the capped distance query and tracer radius.
            keep=(np.maximum(scene[:,:2],scene[:,2:4])+scene[:,4,None]+CAP>=-HALF).all(axis=1)&(np.minimum(scene[:,:2],scene[:,2:4])-scene[:,4,None]-CAP<=HALF).all(axis=1)
            np.savetxt(dest/'solid_roi.tsv',scene[keep],fmt='%.17g')
            region=np.stack(regions[t]);dx=selected[0]['dx_um']
            coords=-23.8+(np.arange(256)+.5)*dx
            # Require every raster cell intersecting the proposed ROI to be inside support.
            roi_cell=(coords+dx/2>=-HALF)&(coords-dx/2<=HALF)
            roi_in_support=bool(region[:,roi_cell][:,:,roi_cell].all())
            np.savez_compressed(dest/'saved_support.npz',region=region,H_um=H,dx_um=dx,
                axis_order='z,y,x',extent_um=[-23.8,23.8,-23.8,23.8],layer_seeds=stack['layer_seeds'])
            mean_area=st['mean_saved_support_area_um2'];vol=50*H*mean_area
            g=dict(state_id=state_id,stack_id=sid,replicate=stack['replicate'],f0_nN=stack['f0_nN'],
                mechanical_time_s=t,layers=50,H_um=H,Lz_um=L,mean_support_area_um2=mean_area,
                analyzed_volume_um3=vol,C_A=st['C_A_from_saved_masks'],
                saved_mask_phi_geom=st['saved_mask_geometric_porosity'],
                saved_mask_void_volume_um3=vol*st['saved_mask_geometric_porosity'],
                local_ROI_inside_all_saved_supports=roi_in_support)
            globalrows.append(g)
            info=dict(**g,layers_source=metadata[t],full_primitives=len(scene),roi_primitives=int(keep.sum()),
                full_scene_sha256=digest(dest/'solid_full.tsv'),roi_scene_sha256=digest(dest/'solid_roi.tsv'))
            save(dest/'state.json',info);states.append(info)
        print('PREPARED',sid,'3 states',flush=True)
    assert len(states)==24 and all(q['local_ROI_inside_all_saved_supports'] for q in globalrows)
    # Cross-force initial scenes must agree; the f0 comparison starts with the same geometry.
    for rep in (1,2):
        s=[q for q in states if q['replicate']==rep and q['mechanical_time_s']==0]
        assert len({q['full_scene_sha256'] for q in s})==1
    csvsave(OUT/'global_geometry.csv',globalrows)
    save(OUT/'states.json',states)
    plan=dict(status='prepared_geometry; transport not launched; domain choice pending',
        mechanical_times_s=[0,5,10],f0_nN=[0,2.2,4.4,8.8],mechanical_replicates=2,layers=50,H_um=H,Lz_um=L,
        tracer_radius_um=RADIUS,tracer_diameter_um=2*RADIUS,D_um2_s=D,
        candidate_transport=dict(domain='central ROI only; whole-clot boundary method not yet configured',
            xy_bounds_um=[-HALF,HALF],z_bounds_um=[0,L],dt_s=.000015625,horizon_s=100,
            reporting_horizons_s=[10,20,50,100],particles_per_state=N_INITIAL,batch_particles=BATCH_N,
            batches_per_state=N_INITIAL//BATCH_N,maximum_workers=4,executor=str(TRANSPORT_EXE),
            method=6,executor_sha256=digest(TRANSPORT_EXE),
            scope='same collision kernel and ROI as prior reference; new states and native initialization extend its tested scope',
            initial='Uniform proposals in xy disk radius1.19um and z=[2,3]um, conditioned ONLY on continuous solid clearance; disconnected or unresolved starts retained.',
            initial_pairing='Same proposal pool across f0 and times within replicate; accepted positions may differ. Exact initial positions coincide across f0 only when geometry does.',
            trajectory_pairing='Same per-batch RNG seeds and particle index across f0/time; retain candidate ID and actual xyz. No claim of identical starts when they differ.',
            stop_rule='Fixed1024 histories/state, no interim significance-based enlargement. This is a proposed pilot allocation, not a power guarantee.',
            boundaries='Centers absorbed at z=0 or L; lateral or solid hit rejects full proposal and consumes dt; obstacle/lateral priority for ties.',
            physical_scope='Extruded backbone capsules and platelet disks; cohesion acts in mechanics but is not a new solid tube, matching original raster and validation.',
            free_reference='Use actual initial z: mean z0*(L-z0)/(2D), not automatically L^2/(8D). Finite-horizon reference must also use the same censoring.'),
        local_geometry=dict(monte_carlo_points=MC_N,random_seed_rule='[916001,replicate]',
            hxy_um=2*HALF/SIDE,hz_um=L/NZ,graph='Lipschitz inner6 and outer26 cell graphs, connectivity interval; not a converged mask',
            inequality='connected lower <= continuous F_conn <= possible upper for the fixed ROI; unresolved != disconnected'),
        statistics=dict(unit='mechanical stack; two paired replicates/f0, not layers or tracers',
            report='Each replicate, mean and descriptive sample SD; paired changes from0s and f0=0 within replicate.',
            survival='S(t)=mean((event==0) or observed_time>t), t<=Tc; RMST(T)=mean(min(observed_time,T)); censor at100.',
            uncertainty='Tracer and geometry-MC uncertainty conditional on stack; do not present it as between-stack confidence.',
            connectivity='Store certified_connected/certified_disconnected_within_ROI/unresolved per initial center; no exact mixture decomposition while unresolved starts remain.',
            inference='No significance or biological-validation claim with two stacks; no permanent-capture or Ragg inference.'),
        old_transport_pooled=False,mechanics_rerun=False,manuscript_changed=False)
    save(OUT/'protocol.json',plan)
    assert all(digest(p)==h for p,h in protected.items()),'Source mutation'
    save(OUT/'preparation_manifest.json',dict(completed=True,source_sha256=protected,
        protected_inputs_unchanged=True,script_sha256=digest(__file__),states=24,
        initial_scene_pairing_passed=True,elapsed_s=time.perf_counter()-tick))


def measure_local():
    assert (OUT/'preparation_manifest.json').exists()
    assert not (OUT/'local_geometry_manifest.json').exists(),'Completed: no overwrite'
    tests();tick=time.perf_counter();states=json.loads((OUT/'states.json').read_text()); rows=[];jobs=[]
    shared=OUT/'shared';shared.mkdir(exist_ok=True)
    for rep in (1,2):
        rng=np.random.default_rng([916001,rep])
        xyz=rng.uniform([-HALF,-HALF,0],[HALF,HALF,L],(MC_N,3))
        np.savetxt(shared/f'volume_rep{rep}.tsv',xyz,fmt='%.17g')
        rng=np.random.default_rng([916002,rep]);ang=rng.uniform(0,2*np.pi,CANDIDATE_N)
        rho=1.19*np.sqrt(rng.random(CANDIDATE_N));z=rng.uniform(2,3,CANDIDATE_N)
        xyz=np.c_[rho*np.cos(ang),rho*np.sin(ang),z]
        np.savetxt(shared/f'release_candidates_rep{rep}.tsv',xyz,fmt='%.17g')
    for st in states:
        sid=st['state_id'];rep=st['replicate'];dest=OUT/'states'/sid
        assert digest(dest/'solid_roi.tsv')==st['roi_scene_sha256']
        dist=points(dest/'solid_roi.tsv',shared/f'volume_rep{rep}.tsv',dest/'volume_clearance.tsv')
        assert len(dist)==MC_N and np.all((dist>=0)&(dist<=CAP))
        void=dist>EPS;acc=dist>=RADIUS+EPS;assert np.all(~acc|void)
        dstart=points(dest/'solid_roi.tsv',shared/f'release_candidates_rep{rep}.tsv',dest/'candidate_clearance.tsv')
        ids=np.flatnonzero(dstart>=RADIUS+EPS)[:N_INITIAL];assert len(ids)==N_INITIAL
        candidates=np.loadtxt(shared/f'release_candidates_rep{rep}.tsv'); xyz=candidates[ids]
        shell([EXE,'grid',dest/'solid_roi.tsv',dest/'cell_clearance.bin',SIDE,NZ,HALF,L,CAP])
        field=np.fromfile(dest/'cell_clearance.bin',dtype='<f8').reshape(NZ,SIDE,SIDE)
        status,conn=classification(field,xyz)
        row=dict(state_id=sid,replicate=rep,f0_nN=st['f0_nN'],mechanical_time_s=st['mechanical_time_s'],
            roi_volume_um3=(2*HALF)**2*L,phi_geom_MC=float(void.mean()),phi_acc_MC=float(acc.mean()),
            void_volume_um3=float(void.mean()*(2*HALF)**2*L),accessible_volume_um3=float(acc.mean()*(2*HALF)**2*L),
            phi_acc_MC95_low=wilson(int(acc.sum()),MC_N)[0],phi_acc_MC95_high=wilson(int(acc.sum()),MC_N)[1],
            sampled_points=MC_N,initial_count=N_INITIAL,initial_candidate_acceptance=float(np.mean(dstart>=RADIUS+EPS)),
            free_unrestricted_mean_exit_time_s=float(np.mean(xyz[:,2]*(L-xyz[:,2])/(2*D))),**conn)
        rows.append(row)
        np.savez_compressed(dest/'local_geometry_samples.npz',volume_clearance_um=dist,
            initial_xyz_um=xyz,initial_candidate_id=ids,initial_clearance_um=dstart[ids],connectivity_status=status,
            status_definitions='1=connected,-1=disconnected_within_ROI,0=unresolved',
            H_um=H,roi_bounds_xyz_um=[[-HALF,HALF],[-HALF,HALF],[0,L]])
        for b in range(N_INITIAL//BATCH_N):
            a=b*BATCH_N; initial=dest/f'initial_batch{b:02d}.tsv'
            np.savetxt(initial,xyz[a:a+BATCH_N],fmt='%.17g')
            output=dest/f'transport_batch{b:02d}'
            job=dict(state_id=sid,replicate=rep,f0_nN=st['f0_nN'],mechanical_time_s=st['mechanical_time_s'],
                batch=b,particles=BATCH_N,seed=916100+(rep-1)*100+b,status='not_run',
                initial_path=str(initial),initial_sha256=digest(initial),scene_sha256=digest(dest/'solid_roi.tsv'),
                command=[str(TRANSPORT_EXE),'run',str(dest/'solid_roi.tsv'),str(initial),str(output),'100',str(BATCH_N),str(916100+(rep-1)*100+b),str(RADIUS),str(D)])
            jobs.append(job)
        save(dest/'local_metrics.json',row)
        csvsave(OUT/'local_geometry.csv',rows)
        print('GEOMETRY',sid,'acc',round(row['phi_acc_MC'],4),'F_conn bounds',round(conn['F_conn_lower'],3),round(conn['F_conn_upper'],3),flush=True)
    # Identical initial state must yield identical geometry and starts across f0.
    for rep in (1,2):
        s=[r for r in rows if r['replicate']==rep and r['mechanical_time_s']==0]
        for key in ['phi_geom_MC','phi_acc_MC','F_conn_lower','F_conn_upper']:
            assert len({r[key] for r in s})==1
    manifest=json.loads((OUT/'preparation_manifest.json').read_text())
    assert all(digest(p)==h for p,h in manifest['source_sha256'].items())
    unique={}
    for job in jobs:
        key=(job['scene_sha256'],job['initial_sha256'],job['seed'])
        job['job_id']=f"{job['state_id']}_batch{job['batch']:02d}"
        if key in unique:
            job['reuse_job_id']=unique[key]
            job['status']='reuse_planned_identical_geometry_starts_noise'
        else:
            unique[key]=job['job_id']
            job['reuse_job_id']=None
    assert len(unique)==72, 'Expected 18 unique geometries x4 batches'
    save(OUT/'transport_jobs_PROPOSED_NOT_RUN.json',jobs)
    save(OUT/'local_geometry_manifest.json',dict(completed=True,states=len(rows),
        proposed_state_batches=len(jobs),unique_transport_jobs=len(unique),
        redundant_state_batches_reused=len(jobs)-len(unique),
        unique_proposed_histories=len(unique)*BATCH_N,transport_histories_generated=0,mechanics_rerun=False,
        protected_inputs_unchanged=True,executor_sha256=digest(EXE),script_sha256=digest(__file__),
        elapsed_s=time.perf_counter()-tick,domain='local candidate; selection pending'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','local']);a=p.parse_args()
    {'prepare':prepare,'local':measure_local}[a.mode]()
