# Model and scope

Mechanics is simulated independently in 2D layers. Their ordered extrusions
form the transport stack; layers do not exchange mechanical forces. Backbone
and cohesive segments both contribute to the current solid obstacles.

The paper's frozen transport protocol uses 50 layers, H = 0.1 micrometre,
Lz = 5 micrometres, tracer radius 0.25 micrometre, D = 0.4 micrometre²/s,
dt = 0.000015625 s and observation horizon 100 s. The release region is central,
conditioned on geometric validity, not on connection to an exit. Both axial
faces absorb tracer centres. Solid and common lateral boundaries reflect.

The continuous-solid distance test is distinct from the conservative auxiliary
grid used to certify connectivity. An unresolved connection is not a disconnected
tracer. A numerical failure is not censoring.

The current production cohort has nine mechanical seed groups and 117 unique
geometries (144 nominal force/time labels). Histories from the shared initial
states must not be counted four times. Bootstrap intervals are conditional on
the geometry; min-max ranges across networks are not confidence intervals.

The larger temporal study compares three frozen states from two mechanical
seed groups, with 8192 histories per state and step, at three step sizes.
All primary intervals fit ±2.5 pp for S(100 s) and ±3% for relative RMST.
The survival tolerance was revised after inspecting results; the original ±2 pp
criterion remains inconclusive in one state. This is not proof of a zero-step
limit, all-state equivalence, evolving-network validity, or biological calibration.

Effective drags, planar platelet density/stereology, and several adaptation
parameters remain scientific questions. Do not silently mix future revised
mechanics with the existing paper cohort.
