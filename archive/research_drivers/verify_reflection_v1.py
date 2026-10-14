"""Bounded verification of isolated reflection code, with preserved source data."""
from pathlib import Path
import argparse,csv,hashlib,json,subprocess,time
import numpy as np
from scipy.optimize import minimize_scalar,brentq

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/verificacion_reflexion_20260919'
CPP=ROOT/'work/transport_reflection_v1.cpp'
EXE=OUT/'transport_reflection_v1.exe'
SOURCE=ROOT/'outputs/analisis_homogeneo_con_cohesion_20260917'
DOMAINS=ROOT/'outputs/dominio_comun_y_bordes_20260917'

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def call(args):return subprocess.run(list(map(str,args)),check=True,capture_output=True,text=True)
def array(p,x):np.savetxt(p,x,fmt='%.17g')
def domain(path,mask,xmin,ymin,dx,L):
    with path.open('w') as f:
        f.write(f'{mask.shape[1]} {mask.shape[0]} {xmin:.17g} {ymin:.17g} {dx:.17g} {L:.17g}\n')
        np.savetxt(f,mask.astype(int),fmt='%d')
def runrows(mode,scene,reg,rows,name):
    q=OUT/'controls'/f'{name}_input.tsv';o=q.with_name(f'{name}_output.tsv');array(q,rows)
    args=[EXE,mode]+([scene,reg] if mode=='steps' else [reg] if mode=='lateral' else [scene])+[q,o]
    call(args);return np.loadtxt(o,ndmin=2)
def analytic_ball(a,d,center,R,xy=False):
    a=a-center
    if xy:a=a[:2];d=d[:2]
    aa=d@d;bb=2*a@d;cc=a@a-R*R;disc=bb*bb-4*aa*cc
    if disc<0 or aa==0:return 2.
    t=(-bb-np.sqrt(disc))/(2*aa)
    return t if 0<=t<=1 else 2.
def capsule_distance(row,p):
    a=row[:2];v=row[2:4]-a;vv=v@v
    u=np.clip((p[:2]-a)@v/vv,0,1) if vv else 0.
    dxy=max(np.linalg.norm(p[:2]-a-u*v)-row[4],0.)
    dz=max(row[5]-p[2],p[2]-row[6],0.)
    return np.hypot(dxy,dz)
def check():
    OUT.mkdir(exist_ok=True);(OUT/'controls').mkdir(exist_ok=True)
    call(['C:/Strawberry/c/bin/g++.exe','-std=c++17','-O3','-DNDEBUG','-static',CPP,'-o',EXE])
    rng=np.random.default_rng(919001);results=[];tests=OUT/'controls'
    empty=tests/'empty.tsv';array(empty,np.empty((0,7)))
    broad=tests/'broad.txt';domain(broad,np.ones((50,50),bool),-25,-25,1,5)
    for name,primitive,R,xy in [('sphere',[0,0,0,0,0,2.5,2.5],.25,False),
                              ('cylinder',[0,0,0,0,.5,-10,10],.75,True)]:
        scene=tests/f'{name}.tsv';array(scene,[primitive]);rows=[];truth=[]
        center=np.array([0,0,2.5])
        for i in range(1500):
            direction=rng.normal(size=3)
            if xy:direction[2]=0
            a=center+direction/np.linalg.norm(direction)*rng.uniform(R+.0001,2)
            d=rng.normal(size=3)*2
            rows.append(np.r_[a,d,.25]);truth.append(analytic_ball(a,d,center,R,xy))
        got=runrows('contacts',scene,None,rows,name)
        assert np.array_equal(got[:,0]<=1,np.array(truth)<=1),name
        assert np.max(np.abs(got[:,0]-truth))<2e-8,name
        results.append(dict(test=name+'_analytic_first_contact',cases=len(rows),passed=True,
                            max_fraction_error=float(np.max(np.abs(got[:,0]-truth)))))
    primitive=np.array([-.7,-.3,.8,.5,.12,1.1,2.2]);scene=tests/'finite_capsule.tsv';array(scene,[primitive])
    rows=[];truth=[]
    while len(rows)<500:
        a=rng.uniform([-2,-2,.2],[2,2,4.8]);d=rng.normal(size=3)*1.5
        if capsule_distance(primitive,a)<=.25:continue
        f=lambda t:capsule_distance(primitive,a+t*d)-.25
        opt=minimize_scalar(f,bounds=(0,1),method='bounded',options={'xatol':1e-13})
        tm=min([0.,1.,float(opt.x)],key=f)
        t=brentq(f,0,tm,xtol=1e-14) if f(tm)<0 else 2.
        rows.append(np.r_[a,d,.25]);truth.append(t)
    got=runrows('contacts',scene,None,rows,'finite_capsule')
    assert np.array_equal(got[:,0]<=1,np.array(truth)<=1)
    assert np.max(np.abs(got[:,0]-truth))<2e-8
    results.append(dict(test='finite_extruded_capsule_independent_optimizer',cases=500,passed=True,
                        max_fraction_error=float(np.max(np.abs(got[:,0]-truth)))))
    # Planar solid top face: exactly solvable normal reflection, no lateral hit.
    slab=tests/'slab.tsv';array(slab,[[0,0,0,0,100,0,1]])
    noise=rng.normal(0,.1,(20000,3));a=np.tile([0,0,1.26],(len(noise),1))
    got=runrows('steps',slab,broad,np.c_[a,noise,np.full(len(a),.25)],'plane_distribution')
    exact=a+noise;exact[:,2]=1.25+np.abs(exact[:,2]-1.25)
    assert np.max(np.abs(got[:,:3]-exact))<2e-11
    assert np.all(got[:,3]==0)
    results.append(dict(test='planar_reflection_exact_image_per_sample',cases=len(a),passed=True,
                        max_coordinate_error_um=float(np.max(np.abs(got[:,:3]-exact)))))
    # Rectangle exact fold, multiple lateral collisions including a corner.
    rect=tests/'rectangle.txt';domain(rect,np.ones((2,2),bool),-1,-1,1,5)
    rows=np.array([[0,0,2.5,7.3,-5.1,0,.25],[0,0,2.5,3.2,3.2,0,.25],
        [.5,.5,2.5,1,0,0,.25],[0,0,2.5,0,0,4,.25]])
    got=runrows('steps',empty,rect,rows,'rectangle_folds')
    expected=rows[:3,:3]+rows[:3,3:6]
    expected[:,:2]=1-np.abs(np.mod(expected[:,:2]+1,4)-2)
    assert np.max(np.abs(got[:3,:3]-expected))<2e-11,(got,expected)
    assert got[3,3]==2 and abs(got[3,6]-.625)<1e-12
    results.append(dict(test='rectangular_multireflection_corner_and_exit',cases=4,passed=True))
    # Z reflection before top/bottom exit, and obstacle lying beyond first exit.
    rows=np.array([[0,0,1.5,.3,.4,-1,.25],[0,0,1.5,0,0,-5,.25]])
    got=runrows('steps',slab,broad,rows,'solid_exit_order')
    assert np.allclose(got[0,:3],[.3,.4,2.],atol=1e-11) and got[0,3]==0
    assert got[1,3]==2 and got[1,4]>=1 # first reflect from slab, then leave top
    beyond=tests/'beyond.tsv';array(beyond,[[0,0,0,0,100,6,7]])
    got2=runrows('steps',beyond,broad,[[0,0,2.5,0,0,8,.25]],'exit_before_solid')
    assert got2[0,3]==2 and got2[0,4]==0
    results.append(dict(test='solid_before_exit_and_exit_before_solid',cases=3,passed=True))
    # Independent segment-box reference for a nonconvex lateral mask with a hole.
    mask=np.ones((8,8),bool);mask[2:4,3:5]=False;mask[5:,6:]=False
    irregular=tests/'irregular.txt';domain(irregular,mask,-4,-4,1,5)
    cells=np.argwhere(mask);rows=[];truth=[]
    for _ in range(2000):
        y,x=cells[rng.integers(len(cells))];a=np.array([x-4,y-4,2.5])+np.r_[rng.uniform(.001,.999,2),0]
        d=np.r_[rng.uniform(-10,10,2),0];t=2.
        for j in range(2):
            end=a[j]+d[j]
            if end<-4:t=min(t,(-4-a[j])/d[j])
            if end>4:t=min(t,(4-a[j])/d[j])
        for yy,xx in np.argwhere(~mask):
            low=np.array([xx-4,yy-4]);high=low+1
            t1=(low-a[:2])/d[:2];t2=(high-a[:2])/d[:2]
            enter=max(0.,np.minimum(t1,t2).max());leave=min(1.,np.maximum(t1,t2).min())
            if enter<=leave:t=min(t,enter)
        rows.append(np.r_[a,d]);truth.append(t)
    got=runrows('lateral',None,irregular,rows,'lateral_box_reference')
    assert np.max(np.abs(got[:,0]-truth))<2e-11
    results.append(dict(test='lateral_DDA_vs_independent_segment_box',cases=2000,passed=True))
    # Reversibility checks full piecewise paths on curved obstacle + outer walls.
    arena=tests/'arena.txt';domain(arena,np.ones((6,6),bool),-3,-3,1,5)
    cylinder=tests/'cylinder.tsv';rows=[]
    while len(rows)<1000:
        a=np.r_[rng.uniform(-2.9,2.9,2),2.5]
        if np.linalg.norm(a[:2])<.751:continue
        d=np.r_[rng.normal(size=2)*3,0];rows.append(np.r_[a,d,.25])
    rows=np.asarray(rows);got=runrows('steps',cylinder,arena,rows,'forward_curved')
    back=runrows('steps',cylinder,arena,np.c_[got[:,:3],-got[:,7:10],np.full(len(got),.25)],'backward_curved')
    err=np.max(np.abs(back[:,:3]-rows[:,:3]))
    assert err<1e-8,err
    assert np.min(got[:,11])>=.25-2e-10
    results.append(dict(test='curved_and_lateral_multihit_reversibility',cases=1000,passed=True,max_error_um=float(err)))
    # Export selected real domain without changing native support resolution.
    state='rep1_f0_2.2_tm10';source=DOMAINS/'states'/f'{state}.npz'
    with np.load(source,allow_pickle=False) as data:
        domain(OUT/'real_domain.txt',data['common_support'],float(data['extent_um'][0]),
            float(data['extent_um'][2]),float(data['dx_um']),float(data['Lz_um']))
    save(OUT/'controls.json',dict(passed=True,tests=results,total_cases=sum(v['cases'] for v in results),
        source_sha256={str(CPP):sha(CPP),str(ROOT/'work/transport_geometry_bridge.cpp'):sha(ROOT/'work/transport_geometry_bridge.cpp')},
        domain_source=str(source),domain_sha256=sha(source),exe_sha256=sha(EXE)))
    print(json.dumps(results,indent=2))

def smoke():
    state=SOURCE/'states/rep1_f0_2.2_tm10';scene=state/'solid_full.tsv';initial=state/'initial_batch00.tsv'
    folder=OUT/'smoke';folder.mkdir(exist_ok=True)
    for method in ['reflect','reject']:
        started=time.perf_counter()
        res=call([EXE,'run',scene,OUT/'real_domain.txt',initial,folder/f'{method}.csv',8,1,.00003125,919002,.25,.4,int(method=='reject')])
        (folder/f'{method}.log').write_text(res.stdout,encoding='utf-8');print(method,time.perf_counter()-started,res.stdout,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['controls','smoke']);args=parser.parse_args()
    {'controls':check,'smoke':smoke}[args.mode]()
