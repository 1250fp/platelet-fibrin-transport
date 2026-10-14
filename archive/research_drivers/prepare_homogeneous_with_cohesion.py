"""Correct the prepared solid geometry using saved time-dependent cohesion.

No mechanical integration, transport execution, production-code changes or
overwrites of the previous analysis. Geometry is regenerated in a new directory.
"""
from __future__ import annotations
import argparse
import importlib
import json
import sys
import time
import types
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import prepare_homogeneous_analysis as base

ROOT=base.ROOT
OLD=ROOT/'outputs/analisis_homogeneo_20260916'
OUT=ROOT/'outputs/analisis_homogeneo_con_cohesion_20260917'
CORE=Path('E:/tesisprog/fibrina/sub_fibrin_model/program19_3/program19_core')
base.OUT=OUT
# Do not rebuild or replace any binary used by the preceding analysis.
base.EXE=OUT/'homogeneous_clearance.exe'


def layers_module():
    name='_cohesion_geometry_readonly'
    package=types.ModuleType(name);package.__path__=[str(CORE)];sys.modules[name]=package
    return importlib.import_module(name+'.layers')


def capsule_rows(nodes,pairs,r,k):
    pairs=np.asarray(pairs,dtype=int).reshape(-1,2)
    a,b=nodes[pairs[:,0]],nodes[pairs[:,1]]
    return np.c_[a,b,np.full(len(a),r),np.full(len(a),k*base.H),np.full(len(a),(k+1)*base.H)]


def prepare():
    assert not OUT.exists(),'Output exists: inspect it; do not overwrite'
    OUT.mkdir(parents=True);tick=time.perf_counter()
    old_states=base.read(OLD/'states.json');old_prep=base.read(OLD/'preparation_manifest.json')
    assert all(base.digest(p)==h for p,h in old_prep['source_sha256'].items())
    for p in (Path(__file__),Path(base.__file__),base.CPP,base.OLD_CPP,base.TRANSPORT_EXE):base.protect(p)
    for name in ('layers.py','config.py','state.py'):base.protect(CORE/name)
    raster=layers_module().rasterize_geometry
    grouped=defaultdict(list)
    for st in old_states:grouped[st['stack_id']].append(st)
    states=[];global_rows=[];link_rows=[];n_empty_checked=0
    for sid, old_group in grouped.items():
        old_group=sorted(old_group,key=lambda q:q['mechanical_time_s'])
        additions={q['mechanical_time_s']:[] for q in old_group}
        solids={t:[] for t in additions};regions={t:[] for t in additions};meta={t:[] for t in additions}
        for li in range(50):
            source=old_group[0]['layers_source'][li];path=Path(source['path'])
            cfg=base.read(path/'config.json');config=SimpleNamespace(**cfg)
            with np.load(base.protect(path/'trajectory.npz'),allow_pickle=False) as tr:
                nodes,platelets,backbone=tr['nodes_um'],tr['platelets_um'],tr['backbone_pairs']
                offsets,pairs,birth=tr['cohesion_frame_offsets'],tr['cohesion_pairs_flat'],tr['cohesion_birth_times_s_flat']
                with np.load(base.protect(path/'intermediate_masks.npz'),allow_pickle=False) as mh:
                    old_solids,old_regions=mh['solid_masks'],mh['analyzed_region_masks']
                for old_st in old_group:
                    t=old_st['mechanical_time_s'];entry=old_st['layers_source'][li];j=entry['frame_index']
                    cp=pairs[offsets[j]:offsets[j+1]]
                    assert len(cp)==entry['cohesive_links_saved']
                    assert np.all(birth[offsets[j]:offsets[j+1]]<=t+1e-8),'Link before its birth'
                    canonical={tuple(sorted(map(int,e))) for e in cp}
                    assert len(canonical)==len(cp) and not canonical.intersection(tuple(sorted(map(int,e))) for e in backbone)
                    combined=np.vstack([backbone,cp])
                    solid,region,extent,dx=raster(nodes[j],platelets[j],combined,config)
                    assert not np.any(old_solids[j]&~solid),'Adding links removed solid'
                    assert not np.any(old_regions[j]&~region),'Rebuilt support lost original cells'
                    if len(cp)==0:
                        assert np.array_equal(old_solids[j],solid) and np.array_equal(old_regions[j],region)
                        n_empty_checked+=1
                    additions[t].append(capsule_rows(nodes[j],cp,cfg['fiber_radius_um'],li))
                    solids[t].append(solid);regions[t].append(region)
                    m=dict(entry,cohesion_is_geometric_obstacle=True,cohesive_radius_um=cfg['fiber_radius_um'],
                        cohesive_primitives=len(cp),solid_pixels_added=int((solid&~old_solids[j]).sum()),
                        support_pixels_added=int((region&~old_regions[j]).sum()))
                    meta[t].append(m)
                    link_rows.append(dict(stack_id=sid,mechanical_time_s=t,layer=li+1,seed=entry['seed'],
                        cohesive_segments=len(cp),radius_um=cfg['fiber_radius_um'],
                        solid_pixels_added=m['solid_pixels_added'],support_pixels_added=m['support_pixels_added']))
        initial_area=float(np.stack(regions[0.]).sum(axis=(1,2)).mean()*dx*dx)
        for old_st in old_group:
            t=old_st['mechanical_time_s'];state_id=old_st['state_id'];dest=OUT/'states'/state_id;dest.mkdir(parents=True)
            old_scene_path=base.protect(OLD/'states'/state_id/'solid_full.tsv')
            old_scene=np.loadtxt(old_scene_path);extra=np.concatenate(additions[t]);scene=np.vstack([old_scene,extra])
            np.savetxt(dest/'solid_full.tsv',scene,fmt='%.17g')
            np.savetxt(dest/'cohesive_segments.tsv',extra,fmt='%.17g')
            keep=(np.maximum(scene[:,:2],scene[:,2:4])+scene[:,4,None]+base.CAP>=-base.HALF).all(axis=1)&(np.minimum(scene[:,:2],scene[:,2:4])-scene[:,4,None]-base.CAP<=base.HALF).all(axis=1)
            np.savetxt(dest/'solid_roi.tsv',scene[keep],fmt='%.17g')
            solid=np.stack(solids[t]);region=np.stack(regions[t]);void=region&~solid
            coords=extent[0]+(np.arange(region.shape[1])+.5)*dx
            cell=(coords+dx/2>=-base.HALF)&(coords-dx/2<=base.HALF)
            assert region[:,cell][:,:,cell].all()
            mean_area=float(region.sum(axis=(1,2)).mean()*dx*dx)
            volume=float(region.sum()*dx*dx*base.H)
            np.savez_compressed(dest/'raster_with_cohesion.npz',solid_layers=solid,analyzed_region_masks=region,
                dx_um=dx,H_um=base.H,extent_um=extent,axis_order='z,y,x',cohesion_included=True)
            g=dict(state_id=state_id,stack_id=sid,replicate=old_st['replicate'],f0_nN=old_st['f0_nN'],
                mechanical_time_s=t,mean_support_area_um2=mean_area,analyzed_volume_um3=volume,
                C_A=1-mean_area/initial_area,phi_geom_raster_with_cohesion=float(void.sum()/region.sum()),
                void_volume_raster_with_cohesion_um3=float(void.sum()*dx*dx*base.H),
                solid_pixels_added=sum(m['solid_pixels_added'] for m in meta[t]),
                support_pixels_added=sum(m['support_pixels_added'] for m in meta[t]),
                cohesive_segments=len(extra),H_um=base.H,Lz_um=base.L,layers=50)
            global_rows.append(g)
            st=dict(**g,layers_source=meta[t],full_primitives=len(scene),roi_primitives=int(keep.sum()),
                full_scene_sha256=base.digest(dest/'solid_full.tsv'),roi_scene_sha256=base.digest(dest/'solid_roi.tsv'))
            base.save(dest/'state.json',st);states.append(st)
        print('COHESION INCLUDED',sid,'0/5/10 s',flush=True)
    assert len(states)==24 and n_empty_checked>=400
    for rep in (1,2):
        assert len({s['full_scene_sha256'] for s in states if s['replicate']==rep and s['mechanical_time_s']==0})==1
    base.csvsave(OUT/'global_geometry_with_cohesion.csv',global_rows)
    base.csvsave(OUT/'links_by_layer_and_time.csv',link_rows)
    base.save(OUT/'states.json',states)
    protocol=base.read(OLD/'protocol.json')
    protocol['status']='cohesion geometry corrected; transport not launched; domain pending'
    protocol['candidate_transport']['physical_scope']='Union of original backbone and existing cohesive segments, same fiber radius, plus platelet disks; all extruded within their own mechanical slab.'
    protocol['source_geometry_version']='cohesive_solid_v1'
    protocol['historical_without_cohesion_directory']=str(OLD)
    protocol['global_geometry']='Rerasterized original plus cohesive segments; analyzed per-layer support regenerated by the original padded-convex-hull rule, not intersection over layers.'
    base.save(OUT/'protocol.json',protocol)
    assert all(base.digest(p)==h for p,h in base.protected.items())
    base.save(OUT/'preparation_manifest.json',dict(completed=True,source_sha256=base.protected,
        protected_inputs_unchanged=True,script_sha256=base.digest(__file__),states=24,
        empty_cohesion_replays_passed=n_empty_checked,solid_addition_monotone=True,
        birth_time_checks_passed=True,cohesion_radius_matches_fibrin=True,
        production_code_changed=False,mechanics_rerun=False,transport_run=False,
        elapsed_s=time.perf_counter()-tick))


def compare():
    rows=[]
    for st in json.loads((OUT/'states.json').read_text()):
        sid=st['state_id'];dest=OUT/'states'/sid;old=OLD/'states'/sid
        # Same evaluation points: adding closed solids cannot increase clearance.
        d=np.loadtxt(dest/'volume_clearance.tsv');before=np.loadtxt(old/'volume_clearance.tsv')
        assert np.all(d<=before+1e-12)
        new_good=d>=base.RADIUS+base.EPS;old_good=before>=base.RADIUS+base.EPS
        assert not np.any(new_good&~old_good)
        initial_new=np.loadtxt(dest/'candidate_clearance.tsv');initial_old=np.loadtxt(old/'candidate_clearance.tsv')
        assert np.all(initial_new<=initial_old+1e-12)
        with np.load(old/'local_geometry_samples.npz',allow_pickle=False) as f:old_ids=f['initial_candidate_id']
        row=dict(state_id=sid,replicate=st['replicate'],f0_nN=st['f0_nN'],mechanical_time_s=st['mechanical_time_s'],
            cohesive_segments=st['cohesive_segments'],old_phi_acc=float(old_good.mean()),new_phi_acc=float(new_good.mean()),
            delta_phi_acc_pp=float(100*(new_good.mean()-old_good.mean())),
            old_initial_centers_now_invalid=int((initial_new[old_ids]<base.RADIUS+base.EPS).sum()),
            support_pixels_added=st['support_pixels_added'],solid_pixels_added=st['solid_pixels_added'])
        if st['mechanical_time_s']==0:
            assert st['cohesive_segments']==0 and np.array_equal(d,before)
        rows.append(row)
    base.csvsave(OUT/'effect_of_including_cohesion.csv',rows)
    manifest=json.loads((OUT/'preparation_manifest.json').read_text())
    assert all(base.digest(p)==h for p,h in manifest['source_sha256'].items())
    base.save(OUT/'cohesion_checks.json',dict(passed=True,geometry_states_compared=24,
        same_sampling_points=True,clearance_never_increased=True,no_newly_admitted_centers=True,
        t0_unchanged=True,originals_unchanged=True,production_transport_histories=0))
    print('COHESION GEOMETRY CHECKS PASSED',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','local','compare']);a=p.parse_args()
    {'prepare':prepare,'local':base.measure_local,'compare':compare}[a.mode]()
