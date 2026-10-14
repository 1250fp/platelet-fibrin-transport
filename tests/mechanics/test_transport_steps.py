"""Geometric/censor regressions independent of a full mechanical simulation."""
from __future__ import annotations

import unittest
from unittest.mock import patch
import numpy as np

from program19_core.transport_steps import (
    ACCEPT, ESCAPE_Z_MIN, ESCAPE_Z_MAX, REJECT_SOLID, REJECT_LATERAL,
    trace_proposals, empirical_survival,
)
from program19_core.transport import first_passage_simulation
from program19_core.evolving_transport import StackHistory, time_dependent_first_passage


class RecordingRng:
    def __init__(self, jumps=None):
        self.scales=[]
        self.jumps=list(jumps or [])

    def integers(self, low, high, size):
        return np.zeros(size,dtype=int)

    def normal(self, mean, sigma, size):
        self.scales.append(float(sigma))
        result=np.zeros(size)
        if self.jumps:
            result[:]=self.jumps.pop(0)
        return result


class TraversalTests(unittest.TestCase):
    def setUp(self):
        self.mask=np.ones((6,6,6),bool)

    def event(self,p,q,spacing=(1.,1.,1.)):
        return trace_proposals(self.mask,np.asarray([p]),np.asarray([q]),spacing)[0][0]

    def test_within_cell_and_zero_jump(self):
        self.assertEqual(self.event([2.5]*3,[2.6,2.4,2.5]),ACCEPT)
        self.assertEqual(self.event([2.5]*3,[2.5]*3),ACCEPT)

    def test_axial_solid_tunnelling(self):
        self.mask[3,:,:]=False
        self.assertEqual(self.event([2.5]*3,[4.5,2.5,2.5]),REJECT_SOLID)

    def test_obstacle_before_z_escape(self):
        self.mask[3,:,:]=False
        self.assertEqual(self.event([2.5]*3,[8.,2.5,2.5]),REJECT_SOLID)

    def test_z_faces_and_endpoint_contact(self):
        self.assertEqual(self.event([2.5]*3,[-1.,2.5,2.5]),ESCAPE_Z_MIN)
        self.assertEqual(self.event([2.5]*3,[6.,2.5,2.5]),ESCAPE_Z_MAX)
        self.assertEqual(self.event([2.5]*3,[0.,2.5,2.5]),ESCAPE_Z_MIN)

    def test_lateral_order_and_corner_escape(self):
        self.assertEqual(self.event([2.5]*3,[8.,20.,2.5]),REJECT_LATERAL)
        self.assertEqual(self.event([2.5]*3,[20.,8.,2.5]),ESCAPE_Z_MAX)
        self.assertEqual(self.event([2.5]*3,[8.,8.,2.5]),REJECT_LATERAL)

    def test_edge_supercover(self):
        self.mask[2,2,3]=False
        self.assertEqual(self.event([2.5]*3,[2.5,4.5,4.5]),REJECT_SOLID)

    def test_three_axis_corner_supercover(self):
        self.mask[2,3,3]=False
        self.assertEqual(self.event([2.5]*3,[4.5]*3),REJECT_SOLID)

    def test_static_grid_plane_supercover(self):
        self.mask[2,2,3]=False
        self.assertEqual(self.event([3.,2.5,2.5],[3.,2.5,4.5]),REJECT_SOLID)

    def test_anisotropic_spacing(self):
        spacing=np.array([.07,.2,.2])
        self.mask[3,:,:]=False
        self.assertEqual(self.event(np.array([2.5]*3)*spacing,
                                    np.array([4.5,2.5,2.5])*spacing,spacing),REJECT_SOLID)

    def test_invalid_and_empty(self):
        with self.assertRaises(ValueError):self.event([2.5]*3,[np.nan,1.,1.])
        a,f=trace_proposals(self.mask,np.empty((0,3)),np.empty((0,3)),np.ones(3))
        self.assertEqual(len(a),0)

    def test_randomized_against_independent_box_intersections(self):
        rng=np.random.default_rng(8201)
        mask=rng.random((6,6,6))>.25
        free=np.argwhere(mask)
        starts=free[rng.integers(0,len(free),size=1000)]+rng.uniform(.05,.95,(1000,3))
        ends=starts+rng.normal(0,5.,starts.shape)
        actual,_=trace_proposals(mask,starts,ends,np.ones(3))
        blocked=np.argwhere(~mask)
        for i,(p,q) in enumerate(zip(starts,ends)):
            d=q-p; candidates=[(1.+1e-9,ACCEPT)]
            for axis in range(3):
                if d[axis]:
                    t=((6. if d[axis]>0 else 0.)-p[axis])/d[axis]
                    if 0<=t<=1:
                        code=(ESCAPE_Z_MAX if d[axis]>0 else ESCAPE_Z_MIN) if axis==0 else REJECT_LATERAL
                        candidates.append((t,code))
            for cell in blocked:
                lo=(cell-p)/d;hi=(cell+1-p)/d
                entry=max(0.,float(np.minimum(lo,hi).max()))
                leave=min(1.,float(np.maximum(lo,hi).min()))
                if entry<=leave:candidates.append((entry,REJECT_SOLID))
            tmin=min(t for t,_ in candidates)
            expected=max(code for t,code in candidates if abs(t-tmin)<1e-10)
            self.assertEqual(int(actual[i]),expected,f"segment {i}: {p} -> {q}")


class ReportingTests(unittest.TestCase):
    def test_survival_at_censor(self):
        np.testing.assert_array_equal(empirical_survival([.5,np.nan,1.],[0.,.5,1.],1.),[1.,2/3,1/3])

    def test_nondivisible_last_step_variance(self):
        rng=RecordingRng()
        s,t,S=first_passage_simulation(np.ones((50,8,8),bool),.2,.07,1.,4,.003,.007,rng)
        np.testing.assert_allclose(np.square(rng.scales),[.006,.006,.002],atol=1e-14)
        self.assertEqual(s['censored_count'],4)
        self.assertEqual(S[-1,1],1.)
        self.assertAlmostEqual(s['MET_restricted_s'],.007)

    def test_censor_does_not_clamp_a_late_escape(self):
        # With zero proposals all particles stay; no fictitious endpoint escape.
        s,t,S=first_passage_simulation(np.ones((50,8,8),bool),.2,.07,0.,20,.003,.007,RecordingRng())
        self.assertTrue(np.all(np.isnan(t)))
        self.assertEqual(S[-1,1],1.)

    def test_barrier_and_outcome_api(self):
        mask=np.ones((50,8,8),bool);mask[15]=False;mask[34]=False
        s,t,S,out=first_passage_simulation(mask,.2,.07,1.,128,.001,.5,
                                          np.random.default_rng(81),return_outcomes=True)
        self.assertEqual(s['escaped_count'],0)
        self.assertEqual(S[-1,1],1.)
        self.assertTrue(np.all(out==0))

    def test_evolving_and_frozen_use_shared_rule_and_short_last_step(self):
        solid=np.zeros((2,50,8,8),bool);solid[:,15]=True;solid[:,34]=True
        history=StackHistory(np.array([0.,.1]),solid,np.ones_like(solid),.2,.07,('synthetic',))
        def prescribed(s,r,*args):return r & ~s,np.zeros(s.shape)
        for mode in ('frozen','evolving'):
            # Single huge proposal attempts to cross a barrier and escape.
            with patch('program19_core.evolving_transport.accessible_centers',side_effect=prescribed):
                s,a=time_dependent_first_passage(history,release_time_s=0.,observation_time_s=.007,
                    particle_radius_um=.1,diffusion_um2_s=1.,particle_count=4,dt_s=.003,
                    rng=RecordingRng([[8.,0.,0.],[8.,0.,0.],[8.,0.,0.]]),geometry_mode=mode)
            self.assertEqual(s['escaped_count'],0)
            self.assertEqual(a['survival_probability'][-1],1.)
            np.testing.assert_allclose(a['trajectory_positions_zyx_um'][0],a['trajectory_positions_zyx_um'][-1])


if __name__=='__main__':unittest.main()
