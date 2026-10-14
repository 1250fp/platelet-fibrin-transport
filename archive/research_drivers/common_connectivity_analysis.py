"""Certified inner/outer connectivity of the full common domain."""
from pathlib import Path
import json, subprocess,time
import numpy as np
from scipy.ndimage import label,generate_binary_structure

EPS=1e-12
def call(args):subprocess.run(list(map(str,args)),check=True,capture_output=True,text=True)
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def write_domain(npz,path):
    with np.load(npz,allow_pickle=False) as d:
        mask=d['common_support'];extent=d['extent_um'];dx=float(d['dx_um']);L=float(d['Lz_um'])
    with path.open('w') as f:
        f.write(f'{mask.shape[1]} {mask.shape[0]} {extent[0]:.17g} {extent[2]:.17g} {dx:.17g} {L:.17g}\n')
        np.savetxt(f,mask.astype(int),fmt='%d')

def exits(mask,neighbors):
    lab,n=label(mask,generate_binary_structure(3,neighbors))
    hit=np.zeros(n+1,bool);hit[np.unique(lab[0])]=True;hit[np.unique(lab[-1])]=True;hit[0]=False
    connected=hit[lab];del lab
    return connected,n

def compute(exe,scene,domain,initial,dest,sub=2,nz=100,r=.25):
    dest.mkdir(exist_ok=True)
    if (dest/'summary.json').exists():return json.loads((dest/'summary.json').read_text())
    start=time.perf_counter()
    header=domain.read_text().splitlines()[0].split();nx,ny=map(int,header[:2]);xmin,ymin,dx,L=map(float,header[2:])
    shape=(nz,ny*sub,nx*sub);spacing=np.array([L/nz,dx/sub,dx/sub]);diag=.5*np.linalg.norm(spacing)
    file=dest/'clearance_grid.bin';call([exe,'grid',scene,domain,file,sub,nz,r])
    assert file.stat().st_size==np.prod(shape)*8
    d=np.memmap(file,mode='r',dtype='<f8',shape=shape)
    xyz=np.loadtxt(initial,ndmin=2);ij=np.floor((xyz[:,[2,1,0]]-[0,ymin,xmin])/spacing).astype(int)
    assert ((ij>=0)&(ij<shape)).all()
    possible=(d>=0)&(d+diag>=r-EPS)
    assert possible[tuple(ij.T)].all(),'Valid initial excluded by outer cells'
    outer,nouter=exits(possible,3);upper=outer[tuple(ij.T)];del outer,possible
    safe=(d>=0)&(d-diag>=r+EPS)
    connected,ninner=exits(safe,1);lower=connected[tuple(ij.T)].copy()
    kind=np.where(lower,1,0).astype(np.uint8);target=np.full_like(xyz,np.nan)
    # Independent straight-segment witnesses connect boundary-cell starts to the safe graph.
    # A failed search is unresolved, never proof of disconnection.
    qs=[];targets=[]
    reach=.5;span=np.ceil(reach/spacing).astype(int)
    for i in np.flatnonzero(~lower & upper):
        lo=np.maximum(ij[i]-span,0);hi=np.minimum(ij[i]+span+1,shape)
        cells=np.argwhere(connected[lo[0]:hi[0],lo[1]:hi[1],lo[2]:hi[2]])+lo
        if not len(cells):continue
        qzyx=(cells+.5)*spacing+[0,ymin,xmin];q=qzyx[:,[2,1,0]]
        ds=np.linalg.norm(q-xyz[i],axis=1);order=np.argsort(ds)
        for j in order[ds[order]<=reach][:32]:
            qs.append(np.r_[i,xyz[i],q[j]]);targets.append(q[j])
    if qs:
        np.savetxt(dest/'bridges_input.tsv',qs,fmt='%.17g')
        call([exe,'bridges',scene,domain,dest/'bridges_input.tsv',dest/'bridges_output.tsv',r])
        result=np.loadtxt(dest/'bridges_output.tsv',ndmin=2)
        assert len(result)==len(qs)
        for j,(idx,valid) in enumerate(result):
            i=int(idx)
            if valid and not lower[i]:lower[i]=True;kind[i]=2;target[i]=targets[j]
    assert np.all(~lower|upper)
    status=np.where(lower,1,np.where(~upper,-1,0)).astype(np.int8)
    np.savez_compressed(dest/'initial_connectivity.npz',xyz_um=xyz,status=status,certificate_kind=kind,
        bridge_target_um=target,definition='1=certified connected,-1=certified disconnected in common domain,0=unresolved')
    np.savez_compressed(dest/'safe_exit_cells.npz',packed_connected=np.packbits(connected),shape=np.array(shape),
        packed_safe=np.packbits(safe),packed_possible=np.packbits((d>=0)&(d+diag>=r-EPS)),
        spacing_zyx_um=spacing,origin_zyx_um=[0,ymin,xmin])
    summary=dict(completed=True,N=len(xyz),F_conn_lower=float(lower.mean()),F_conn_upper=float(upper.mean()),
        connected=int(lower.sum()),disconnected=int((~upper).sum()),unresolved=int((upper&~lower).sum()),
        connected_by_safe_cell=int((kind==1).sum()),connected_by_bridge=int((kind==2).sum()),
        shape_zyx=list(map(int,shape)),hxy_um=float(spacing[1]),hz_um=float(spacing[0]),halfdiag_um=float(diag),
        inner_components=int(ninner),outer_components=int(nouter),bridge_queries=len(qs),wall_s=time.perf_counter()-start,
        graph='inner6 whole-cell certificate;outer26 Lipschitz superset; independent continuous segment bridges',
        scope='Bracketing connectivity, not spatial convergence; no trajectories used to classify starts')
    # Retain compact cell certificates; this temporary file was created above.
    # Do not retain ~200 MB per state when masks suffice to reproduce the graph.
    del d,safe,connected
    file.unlink()
    summary['temporary_distance_grid_removed']=True
    save(dest/'summary.json',summary);return summary

def tests(exe,dest):
    dest.mkdir(exist_ok=True)
    dom=dest/'domain.txt'
    with dom.open('w') as f:f.write('8 8 -1 -1 .25 2\n'+' '.join(['1']*64))
    initial=dest/'initial.tsv';np.savetxt(initial,[[0,0,1],[.4,.3,1]],fmt='%.17g')
    empty=dest/'empty.tsv';empty.write_text('')
    a=compute(exe,empty,dom,initial,dest/'empty',sub=2,nz=20)
    assert a['connected']==2 and a['disconnected']==0
    blocked=dest/'blocked.tsv';np.savetxt(blocked,[[0,0,0,0,100,.3,.4],[0,0,0,0,100,1.6,1.7]],fmt='%.17g')
    b=compute(exe,blocked,dom,initial,dest/'blocked',sub=2,nz=20)
    assert b['disconnected']==2 and b['connected']==0
    queries=dest/'q.tsv';np.savetxt(queries,[[0,0,0,1,0,0,0],[1,0,0,1,.2,.2,1]],fmt='%.17g')
    call([exe,'bridges',blocked,dom,queries,dest/'qout.tsv','.25'])
    q=np.loadtxt(dest/'qout.tsv');assert np.array_equal(q[:,1],[0,1])
    save(dest/'tests.json',dict(passed=True,controls=['empty all exit-connected','closed middle disconnected','blocked and free continuous bridges']))
