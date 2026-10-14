"""Version2 retains failed v1 results and resolves long paths by subdivision."""
import sys
from pathlib import Path
import numpy as np
import verify_reflection_v1 as base

base.OUT=base.ROOT/'outputs/verificacion_reflexion_20260919_v2'
base.CPP=base.ROOT/'work/transport_reflection_v2.cpp'
base.EXE=base.OUT/'transport_reflection_v2.exe'
mode=sys.argv[1]
if mode=='controls':
    base.check()
    old=base.ROOT/'outputs/verificacion_reflexion_20260919'
    scene=base.SOURCE/'states/rep1_f0_2.2_tm10/solid_full.tsv'
    rows=np.loadtxt(old/'problem_step.tsv',ndmin=2)
    got=base.runrows('steps',scene,base.OUT/'real_domain.txt',rows,'narrow_gap_regression')
    reference=np.loadtxt(old/'problem_step_diagnostic.tsv',ndmin=2)
    assert got[0,12]>0
    assert got[0,4]==reference[0,4]==713
    err=float(np.max(np.abs(got[0,:3]-reference[0,:3])))
    assert err<1e-9,(got,reference)
    import json
    p=base.OUT/'controls.json';result=json.loads(p.read_text())
    result['tests'].append(dict(test='narrow_gap_subdivision_vs_713_bounce_replay',cases=1,passed=True,
        max_error_um=err,geometric_subdivisions=int(got[0,12])))
    thin=base.OUT/'controls/thin_rectangle.txt';base.domain(thin,np.ones((2,2),bool),-5e-6,-5e-6,5e-6,5)
    trial=np.array([[1.234e-6,0,2.5,.01012345,0,0,.25],
                    [1.234e-6,0,2.5,.01012345,0,3,.25]])
    folded=base.runrows('steps',base.OUT/'controls/empty.tsv',thin,trial,'thin_rectangle_subdivision')
    half=5e-6
    for i,frac in enumerate([1.,2.5/3]):
        xx=trial[i,0]+trial[i,3]*frac
        exact=half-abs((xx+half)%(4*half)-2*half)
        assert abs(folded[i,0]-exact)<1e-10,(folded,exact)
        assert folded[i,12]>0
    assert folded[1,3]==2 and abs(folded[1,6]-2.5/3)<1e-11
    result['tests'].append(dict(test='thin_rectangle_1000_reflections_and_exit_with_subdivision',cases=2,passed=True))
    result['total_cases']+=3;base.save(p,result)
    print(result['tests'][-2:],flush=True)
elif mode=='run':
    import run_reflection_verification as run
    # Preserve the fixed pilot allocation and seed rule; record the correct CPP.
    run.main()
elif mode=='report':
    import report_reflection_verification as report
    report.main()
else:raise ValueError(mode)
