# Temporal Green-circuit audit

The 54-case Green-circuit reduction passed for simultaneous local input
patterns. This audit asks whether the same reduced object generalizes in time.

No new morphology, conductance fit or transport fit is introduced.

## Locked temporal programs

All three compact Gate-20 sites receive one HUMAN excitatory event with the
following delays relative to the common 60 ms reference:

~~~text
synchronous   [ 0,  0,  0] ms
forward_5     [ 0,  5, 10] ms
reverse_5     [10,  5,  0] ms
spread_15     [ 0, 15, 30] ms
~~~

The site order is the existing proximal-to-distal ordering of each compact
branch.

These four programs are locked before the run.

## Reduced object

Reuse exactly the successful Green circuit:

~~~text
released HUMAN AMPA/NMDA law
        x
3x3 local Green matrix G_g
        x
site-to-soma transport T_g
~~~

The AMPA and voltage-independent NMDA conductance templates are still measured
from **one single original-geometry HUMAN event**. A delayed event is represented
only by shifting that same template by an integer number of 0.025 ms samples.

No temporally patterned full-model current, voltage or soma trace is used to
construct the reduced input.

## Geometry panel

For every branch:

~~~text
1.00 x original length
0.80 x held-out length
1.20 x held-out length
~~~

Diameter, topology, normalized synapse x locations, mechanisms and event
multiplicity remain unchanged.

## Cases

~~~text
6 branches
x 3 geometries
x 4 temporal programs
= 72 full-model comparisons
~~~

## Diagnostics

For every case report:

1. reduced current-waveform NRMSE;
2. reduced soma-trace NRMSE;
3. transport-oracle soma NRMSE using actual full-model currents;
4. fixed-point convergence;
5. somatic spike guard.

Also compute an **open-loop synapse attacker**:

~~~text
J_open(t) = synapse_law(V_noinput(t), timed conductances)
V_open    = T_g[J_open]
~~~

This removes local voltage feedback while keeping the same conductance
templates and transport. It is diagnostic; no result is tuned against it.

## Locked rulers

Earn:

~~~text
TEMPORAL_GREEN_CIRCUIT_GENERALIZES_WITHOUT_REFIT
~~~

only if:

1. median transport-oracle soma NRMSE <= 0.01;
2. median reduced soma NRMSE across all 72 cases <= 0.02;
3. each temporal program has median reduced soma NRMSE <= 0.03;
4. median reduced current-waveform NRMSE <= 0.02;
5. the nonlinear fixed point converges in 72 / 72 cases;
6. no actual full-model run crosses the conservative -20 mV soma guard.

The open-loop attacker is reported but is not part of the pass ruler because
the purpose of this assay is temporal portability of the already-earned
nonlinear circuit, not a new proof that NMDA feedback exists.

If the transport oracle passes but the reduced circuit fails:

~~~text
TRANSPORT_TEMPORALLY_VALID_LOCAL_NONLINEAR_REDUCTION_NOT_PORTABLE
~~~

If the transport oracle itself fails:

~~~text
TEMPORAL_GREEN_TRANSPORT_INADEQUATE
~~~

No event delays, time alignment, conductance gains, damping, branches, geometry
scales or analysis windows are changed after seeing the result.

This is an architecture audit, not Gate 25.


## Receipt — temporal portability passes

The locked 72-case audit passed.

~~~text
6 branches x 3 geometries x 4 timing programs = 72

median transport-oracle soma NRMSE       0.0037
median reduced-circuit soma NRMSE        0.0050
median reduced current-waveform NRMSE    0.0026

timing-family medians
  synchronous                            0.0070
  forward_5                              0.0053
  reverse_5                              0.0054
  spread_15                              0.0041

fixed-point convergence                  72 / 72
actual soma spike guard                  0
~~~

Classification:

~~~text
TEMPORAL_GREEN_CIRCUIT_GENERALIZES_WITHOUT_REFIT
~~~

The same single-event conductance template, local Green matrix construction,
HUMAN magnesium-block law and transport construction therefore generalize to
forward, reverse and widely staggered three-site event sequences without any
temporal fit.

### Nonlinear feedback ablation

The diagnostic open-loop attacker keeps the exact same conductance programs and
the exact same geometry-specific transport, but evaluates synaptic current only
at the no-input local voltage:

~~~text
J_open = synapse_law(V_noinput, timed conductances)
V_open = T_g[J_open]
~~~

Result:

~~~text
median open-loop soma NRMSE              0.4141
median reduced/open-loop error ratio     0.0121
reduced circuit beats open loop          72 / 72
~~~

Thus the high accuracy is not explained by passive transport plus a prescribed
conductance waveform. The self-consistent local voltage feedback around the
NMDA conductance is essential in this operating regime.

GitHub Actions: run 33501182766, job 99834644732.
