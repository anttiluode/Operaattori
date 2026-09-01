# Operaattori — measured causal channels

This page is deliberately shorter than the gate history.

The scaffold line asks what a local matrix motion is physically allowed to
change.

## Current channel table

| channel | Gate | intervention | receipt | classification |
|---|---:|---|---:|---|
| **intrinsic geometry** | 22 | 20% stretch of the selected subtree's intrinsic cable metric | **36.34% median passive-transfer change** | metric is causal |
| **pure isometric embedding** | 22 | 35-degree parent-local SE(3) bend, same lengths/radii/topology | **184.107 um motion; 5.999e-14 max transfer change** | `CABLE_MODEL_IGNORES_ISOMETRIC_3D_EMBEDDING` |
| **world-space sampling + passive readout** | 23 | original/bent material samples the same smooth field; total drive normalized | **2.11% median soma change; 19.9% >5%; max 22.3%** | `FIXED_WORLD_FIELD_EFFECT_WEAK` |
| **nonlinear branch compartments** | 20 | compact HUMAN NMDA interaction vs complete branch-alone cross-branch null | **70.19% within; 2.01% cross; 14/15 positive modularity pairs** | `HUMAN_GAMMA_STRENGTHENS_SEMI_INDEPENDENT_COMPARTMENTS` |
| **real extracellular field x isometric embedding x nonlinear readout** | 24 | same intrinsic cell, 1 V/m field, 35-degree re-embedding | **1.0003x proximal bend factor; 0/3 frequencies >5%; exact zero-field/material-lock controls** | `FIELD_COUPLING_PRESENT_BUT_NOT_COMPARTMENT_SELECTIVE` |

These rows should not be merged into one vague word such as "geometry."

### Symmetry audit of the pure-embedding row

Gate 22's single 35-degree example is now backed by a deterministic
implementation audit:

~~~text
6 bifurcations
72 local rigid bends
3 local axes x 4 signed angles

max visual displacement               619.350 um
max cable-length error                1.652e-13 um
max passive-transfer change           9.710e-14

20% metric positive controls
median passive-transfer change        11.84%
controls >1%                          6 / 6
~~~

Classification:

~~~text
KNOWN_REEMBEDDING_SYMMETRY_REPLICATED_ON_REAL_SCAFFOLD
~~~

This is a replication of a known cable-geometry invariance on the Operaattori
implementation, not a new curvature theorem. See
[SYMMETRY.md](SYMMETRY.md) and
[results/SYMMETRY_AUDIT.md](results/SYMMETRY_AUDIT.md).

## 1. Intrinsic geometry

The cable equation depends on physical cable properties: length, radius,
membrane parameters, axial resistance and topology.

Gate 22's positive control changes the local matrix translations inside one
subtree by a fixed 20%, preserving radii and topology.

Receipt:

~~~text
median affected-path passive-transfer difference = 0.3634
fraction affected paths >1%                      = 1.000
~~~

Reproduce:

~~~bash
python experiments/gate22_embedding_vs_metric.py \
  --paths 32 --twist-deg 35 --stretch 1.20
~~~

## 2. Pure isometric embedding

The same Gate 22 rotates one local matrix frame without changing cable lengths.

Receipt:

~~~text
max distal displacement                184.107 um
max cable-length change                1.483e-13 um
max passive-transfer relative change   5.999e-14
~~~

This is a fence:

> a moving XYZ drawing is not automatically a moving intrinsic electrical
> operator.

## 3. World-space sampling with a passive readout

Gate 23 adds an explicitly world-space object.

The intrinsic cable operator is reused unchanged. Original and bent embeddings
sample the same deterministic family of smooth positive fields. Each embedding
is normalized to the same total injected drive.

Receipt:

~~~text
material-locked max difference         0
median fixed-world soma difference     0.0211
conditions >5%                         0.199
conditions >10%                        0.060
maximum                                0.2230
~~~

Reproduce:

~~~bash
python experiments/gate23_spatial_field.py \
  --directions 48 --twist-deg 35
~~~

Important: this is an abstract sampling model, **not** a uniform extracellular
electric-field implementation.

## 4. Nonlinear branch compartments

Gate 20 uses the authors' pinned released cell-1125 NEURON model.

Its cross-branch null is already the sum of the two complete nonlinear
branch-alone responses. Therefore cross-branch interaction is measured *after*
within-branch nonlinearity has already been retained.

Receipt:

~~~text
HUMAN
median within-branch nonlinearity       0.7019
median cross-branch interaction         0.0201
median modularity margin log            0.4582
pairs margin >= log(1.05)              14 / 15

rest-matched gamma=.062 margin          0.0579
HUMAN extra gamma-specific margin       0.3971
~~~

Reproduce after installing NEURON and checking out FCI commit
`75ad8b4d81a7f51bf888b30650c543592340db06`:

~~~bash
python experiments/gate20_compartment_modularity.py \
  --fci-root /tmp/fci \
  --branches 6 --sites 3 --cluster-span-um 55 --multiplicity 8
~~~

## What is potentially special here?

Not the individual ingredients.

Matrix exponentials, cable theory, morphology-dependent field polarization and
electric-field effects on NMDA integration all have prior literature.

The potentially useful object is the **causal decomposition on one movable real
scaffold**:

~~~text
matrix motion
   |
   +-- changes intrinsic metric --------> electrical operator changes strongly
   |
   +-- pure rigid re-embedding ---------> intrinsic cable unchanged
   |
   +-- changes relation to world -------> sampled drive can change
   |
   +-- sampled drive enters nonlinear
       branch compartments -------------> tested: coupling exists, selective
                                           amplification is tiny here
~~~

Gate 24 now fills that last arrow negatively at the preregistered scale: the
field coupling is real, but the isometric bend changes the nonlinear
interaction by only about 0.03% in the large-subtree arm.

## Literature fence for Gate 24

Aspart, Remme & Obermayer (2018) derive the field sensitivity of straight and
bent passive cables. In a uniform field, extracellular potential along the bent
cable follows the cable's projection onto the field axis, and bending can
produce frequency-dependent polarization.

https://doi.org/10.1371/journal.pcbi.1006124

Fan et al. (2023; print 2024) explicitly model weak-electric-field regulation
of NMDA dendritic integration and NMDA-spike generation.

https://doi.org/10.1007/s11571-022-09922-y

Therefore Gate 24 is **not** allowed to claim that field-modulated NMDA spikes
are new.

Its narrower question is whether an **isometric local scaffold motion** on this
real cell, with its already-measured branch compartments, changes compartment
recruitment under one fixed extracellular field while the intrinsic cable is
held constant.

## Rule for future gates

A future gate belongs on the scaffold line only if it:

1. fills a missing causal arrow in the diagram above; or
2. kills one.

Do not add life simulations, growth, agents or learning merely because the
scaffold can move.


## 5. Real extracellular field x nonlinear compartment — Gate 24

Gate 24 inserts a uniform 1 V/m extracellular potential through NEURON rather
than converting a world field into injected current.

~~~text
zero-field control                         0
material-locked control                    0
proximal large-subtree HUMAN bend factor   1.0003x
distal small-subtree HUMAN bend factor     1.0002x
proximal frequencies >5%                   0 / 3
proximal/distal median effect ratio         1.206x
spike guard                                 0
~~~

Classification:

~~~text
FIELD_COUPLING_PRESENT_BUT_NOT_COMPARTMENT_SELECTIVE
~~~

This is useful because it separates two statements that are easy to blur:

1. embedding can become physically relevant once a world-space field exists;
2. that does **not** imply the resulting perturbation is large enough to
   recruit a nonlinear compartment differently.

The locked weak-field assay supports the first and rejects the second.


## Weak-field amplitude audit — route closure

After Gate 24, field amplitude was audited without changing any other knob:

~~~text
field        bend factor
0.25 V/m     1.000069x
0.50 V/m     1.000138x
1.00 V/m     1.000277x
2.00 V/m     1.000559x

through-origin R2                       0.999970
local first-order field scale for 5%    175.049 V/m
~~~

Classification:

~~~text
WEAK_FIELD_BEND_EFFECT_IN_LOCAL_LINEAR_REGIME
~~~

This is not a high-field prediction. It says the measured weak-field effect is
a clean, tiny small-signal coupling. Increasing field strength until the 5%
ruler is crossed is explicitly not part of the protocol.

See [results/FIELD_SCALING_AUDIT.md](results/FIELD_SCALING_AUDIT.md).


## Metric → nonlinear interaction audit

The strong Gate-22 metric channel was connected directly to the six Gate-20
nonlinear compact branches by increasing only each selected section's
intrinsic cable length by 20%.

~~~text
median passive site factor                 1.0877x
passive branches >1%                       6 / 6

median exact nonlinear interaction factor  1.0274x
nonlinear branches >5%                     1 / 6
spike guard                                0
~~~

Classification:

~~~text
METRIC_CHANGES_PASSIVE_TRANSPORT_BUT_NOT_NONLINEAR_RATIO
~~~

The intervention therefore changes transport robustly without comparably
changing the normalized three-site nonlinear superposition law. This keeps
`intrinsic transport` and `local nonlinear law` as distinct causal
channels rather than merging them into one morphology effect.

See
[results/METRIC_NONLINEARITY_AUDIT.md](results/METRIC_NONLINEARITY_AUDIT.md).


## Operator composition — transport x local nonlinearity

The causal channels were finally composed rather than merely compared.

For each of the six compact branches, the original geometry supplies a local
nonlinear operator output: the three site-wise AMPA+NMDA current waveforms.

A separate site-to-soma impulse-response operator is then measured for each
geometry.

Held-out branch-length perturbations:

~~~text
scale = 0.80 or 1.20
12 branch x geometry cases

original reconstruction NRMSE       0.0046
frozen-soma attacker NRMSE          0.0982
factorized NRMSE                    0.0512
transport-oracle NRMSE              0.0046
factorized / frozen                 0.5213
factorized wins                     9 / 12
~~~

Classification:

~~~text
TRANSPORT_X_LOCAL_NONLINEAR_OPERATOR_FACTORIZATION
~~~

The oracle demonstrates that site-current to soma transport is almost perfectly
linear at this operating regime. The remaining error comes mainly from the
fact that changing geometry changes the nonlinear synaptic current waveform
slightly, while the reusable architecture intentionally freezes the original
local operator.

See [ARCHITECTURE.md](ARCHITECTURE.md) and
[results/OPERATOR_FACTORIZATION.md](results/OPERATOR_FACTORIZATION.md).


## Cross-input transport reuse

One T_g per held-out geometry was reused across three different local current
operators on every compact branch:

~~~text
middle single
outer pair
triple
~~~

Across 36 held-out branch x geometry x pattern cases:

~~~text
frozen-soma attacker NRMSE       0.0931
factorized NRMSE                 0.0282
transport-oracle NRMSE           0.0039
factorized / frozen              0.3027
factorized wins                  28 / 36

pattern medians
  middle single                  0.0125
  outer pair                     0.0373
  triple                         0.0512
~~~

Classification:

~~~text
TRANSPORT_OPERATOR_REUSES_ACROSS_INPUT_PATTERNS
~~~

Thus the measured transport module is not tied to the triple-input waveform.
Its error remains dominated by portability of the local nonlinear current
operator as interaction complexity grows.

See
[results/CROSS_INPUT_TRANSPORT_AUDIT.md](results/CROSS_INPUT_TRANSPORT_AUDIT.md).


## Local Green matrix x synapse-law reduction

The per-pattern nonlinear-current lookup was removed.

Reduced object:

~~~text
released HUMAN AMPA/NMDA law
        x
3x3 local current-to-voltage Green matrix
        x
3 site-to-soma transport kernels
~~~

Across 54 branch x geometry x pattern cases:

~~~text
transport oracle soma NRMSE       0.0040
reduced soma NRMSE                0.0043
reduced current NRMSE             0.0038

held-out frozen-current NRMSE     0.0282
held-out reduced NRMSE            0.0043
reduced/frozen                    0.1527
reduced wins                      32 / 36
~~~

Classification:

~~~text
LOCAL_GREEN_MATRIX_X_SYNAPSE_LAW_REDUCES_RELEASED_NEURON
~~~

Thus the current decomposition is not merely transport plus a stored nonlinear
waveform. The local nonlinear current is regenerated from the released synapse
law and the geometry-specific local Green matrix.

See
[results/GREEN_CIRCUIT_AUDIT.md](results/GREEN_CIRCUIT_AUDIT.md).


## Temporal portability of the reduced circuit

The Green circuit was reused without refit for four three-site timing programs
across six branches and three geometries.

~~~text
transport oracle soma NRMSE       0.0037
reduced soma NRMSE                0.0050
reduced current NRMSE             0.0026

open-loop soma NRMSE              0.4141
reduced/open-loop                 0.0121
reduced beats open-loop           72 / 72
~~~

Classification:

~~~text
TEMPORAL_GREEN_CIRCUIT_GENERALIZES_WITHOUT_REFIT
~~~

The result separates the architecture into an extremely accurate linear
transport module and an essential nonlinear local voltage-feedback loop.

See
[results/TEMPORAL_GREEN_CIRCUIT_AUDIT.md](results/TEMPORAL_GREEN_CIRCUIT_AUDIT.md).


## Cross-cell operator boundary

A morphology-only leave-one-cell-out audit was run across all 24 released FCI
morphologies under one matched-passive electrical regime.

~~~text
144 measured branch operators

morphology predictor joint NRMSE    0.3522
nearest training branch             0.3199
training mean                        0.6686
training-basis PCA oracle            0.0307

held-out cells beating nearest       11 / 24
local-only morphology map            0.3695
local + whole-cell map               0.3522
~~~

Classification:

~~~text
CROSS_CELL_OPERATOR_LOW_DIMENSIONAL_BUT_MORPHOLOGY_MAP_WEAK
~~~

This separates two facts that should not be conflated:

~~~text
operator family across cells         highly compressible
simple morphology -> operator chart  inadequate
~~~

The extreme human L5 2057 morphology is retained as a real stress case, but
removing it after the fact does not rescue the conclusion.

See
[results/CROSS_CELL_OPERATOR_AUDIT.md](results/CROSS_CELL_OPERATOR_AUDIT.md).


## Electrotonic cross-cell chart

A morphology-only Rall/electrotonic chart was tested on the same 24-cell LOCO
panel.

~~~text
electrotonic chart joint NRMSE      0.4070
gross chart                         0.3522
combined chart                      0.3062
PCA oracle                          0.0307
~~~

Classification:

~~~text
OPERATOR_LOW_DIMENSIONAL_ELECTROTONIC_CHART_STILL_INSUFFICIENT
~~~

Electrotonic descriptors are complementary to gross morphology but do not
provide a sufficient cross-cell coordinate system.

The likely missing object is distributed cable-tree loading rather than another
single scalar descriptor.

See
[results/ELECTROTONIC_CHART_AUDIT.md](results/ELECTROTONIC_CHART_AUDIT.md).


## Full morphology graph -> passive operator

The scalar cross-cell charts were replaced by a direct physical construction.

For each of all 24 released FCI morphologies, a hand-built matched-passive
compartment graph used:

~~~text
section topology
segment length / diameter / membrane area
Ra = 150 ohm cm
Cm = 1 uF/cm^2
Rm = 20,000 ohm cm^2
~~~

Zero-capacitance branch junctions were eliminated analytically and the resulting
linear system generated the same three-site local G and soma T operators.

~~~text
144 branch operator packs

joint G/T NRMSE                  0.0021
local G                          0.0013
soma T                           0.0024
median cell error                0.0018
cells <= 0.10                    23 / 24
~~~

Classification:

~~~text
MORPHOLOGY_GRAPH_GENERATES_PASSIVE_OPERATOR
~~~

This resolves the earlier low-dimensional-map failure:

~~~text
small scalar morphology chart    inadequate
full loaded cable graph          sufficient
~~~

The result is classical cable physics used as an architectural boundary, not a
novel cable-theory claim.

See
[results/DIRECT_CABLE_GRAPH_AUDIT.md](results/DIRECT_CABLE_GRAPH_AUDIT.md).


## Cross-cell morphology -> nonlinear response

The full morphology graph was composed with the fixed HUMAN_PROBE
AMPA/NMDA law across 24 morphologies, three branches per cell and four timing
programs.

The first global-waveform fixed-point implementation failed portability:

~~~text
median soma NRMSE                0.0042
global fixed-point convergence   207 / 288
~~~

Classification:

~~~text
CROSS_CELL_GRAPH_TRANSPORT_VALID_NONLINEAR_CLOSURE_NOT_PORTABLE
~~~

A locked causal solver diagnosis changed no scientific parameter. It evolved
the passive graph state one time step at a time and solved only the three
current-bearing site voltages implicitly.

~~~text
288 cases

graph-current oracle soma NRMSE  0.00404
causal soma NRMSE                0.00253
causal local-voltage NRMSE       0.00080
causal current NRMSE             0.01722

median cell soma NRMSE           0.00214
cells <= 0.10                    23 / 24
Newton convergence               288 / 288
max Newton iterations            4
~~~

Classification:

~~~text
CAUSAL_MORPHOLOGY_GRAPH_NONLINEAR_CLOSURE_VALID
~~~

Failure mechanism:

~~~text
GLOBAL_WAVEFORM_PICARD_WAS_THE_NONPORTABLE_COMPONENT
~~~

Thus the cross-cell architecture requires preserving causal passive state; a
single global waveform Picard loop is not a universally robust runtime.

See
[results/CAUSAL_NONLINEAR_GRAPH_DIAGNOSIS.md](results/CAUSAL_NONLINEAR_GRAPH_DIAGNOSIS.md).
