# Operaattori

![cell 1125 reconstructed as a local-matrix scaffold](cell1125_original_scaffold.png)

**Operaattori** is Finnish for **operator**.

## [▶ Live Compiler Lab](https://anttiluode.github.io/Operaattori/)

The public page now puts the surviving architecture in one interaction:

~~~text
real morphology
     |
     +--> pose only ---------> compiler no-op
     |
     +--> intrinsic metric --> rebuild passive operator
                                  |
                                  v
                         local implicit NMDA
                                  |
                                  v
                              soma trace
~~~

The left pane is the archived **real human cell-1125** Bend Fence. A 35° rigid
bend moves the selected dendritic subtree **184.107 µm** while the measured
classical passive transfer changes by only **5.999e-14**. Switching to a 20%
intrinsic stretch changes the measured passive transfer by **36.34% median**.

The right pane is a transparent **reference specimen** of the extracted causal
runtime. It recompiles a passive cable from intrinsic metric, solves the
voltage-dependent three-site NMDA step with Newton iteration, plots the soma
trace, and compares it with a frozen-rest-current attacker. Rigid screen-space
pose is deliberately invisible to that compiler; intrinsic stretch rebuilds
the operator live.

The **HUNT NMDA KNEE** control adapts the old variance-seeking feedback idea into
an adversarial validator: small ±2.5% drive probes scan for the regime with the
largest local soma gain. Its purpose is to seek the nonlinear region that a
median error can hide, not to manufacture a new biological measurement.

**Fence:** the browser circuit is illustrative. It is not a redistributed
cell-1125 kernel pack and it does not re-measure the 24-cell NEURON audit.
Those accuracy claims remain attached to the archived receipts.

### [Open the Bend Fence directly](https://anttiluode.github.io/Operaattori/bend_fence.html)

The Bend Fence keeps the locked Gate-22 measurements separate from the browser
geometry. A separate **72-bend symmetry audit** repeated the pose null across
six bifurcations, all three local axes and positive/negative bend angles: the
largest subtree motion was **619.350 µm**, while the largest passive-transfer
change was only **9.710e-14**. The 20% intrinsic-metric controls changed
transfer by **11.84% median**. See
**[results/SYMMETRY_AUDIT.md](results/SYMMETRY_AUDIT.md)**.

## Morphology to nonlinear response across 24 cells

The passive graph can now be composed with the explicit voltage-dependent
AMPA/NMDA law without electrically calibrating the target cell.

A standardized **HUMAN_PROBE** was placed on three deterministic apical
branches of every one of the 24 released morphologies and driven with four
timing programs:

~~~text
24 cells x 3 branches x 4 timings = 288 full NEURON comparisons
~~~

The target prediction starts from morphology only:

~~~text
morphology
   |
   v
hand-built passive cable graph
   |
   v
causal 3-site nonlinear closure
   |
   v
soma response
~~~

Result:

~~~text
graph-current oracle soma NRMSE     0.40%
causal nonlinear soma NRMSE         0.25%
causal local-voltage NRMSE          0.08%
causal current NRMSE                1.72%

median cell soma NRMSE              0.21%
cells <= 10%                        23 / 24

Newton convergence                  288 / 288
max Newton iterations               4
~~~

Classification:
`CAUSAL_MORPHOLOGY_GRAPH_NONLINEAR_CLOSURE_VALID`.

The first cross-cell implementation is retained as an important failed arm: a
global fixed-point solve of the entire Green-kernel waveform converged in only
**207 / 288** cases even though its median soma error was 0.42%.

The diagnosis therefore earned:

`GLOBAL_WAVEFORM_PICARD_WAS_THE_NONPORTABLE_COMPONENT`.

So the scientific architecture survived, but **causal state matters to the
numerical architecture**.

That passing numerical organization is now extracted from the audit harness as
**[reduced/causal_graph_circuit.py](reduced/causal_graph_circuit.py)**: a
pure-NumPy compiled state-space runtime with an implicit local NMDA solve. Its
unit tests cover the backward-Euler compiler, analytic current derivative,
zero-input rest state, nonlinear site closure, and invariance to compartment
relabeling.

The compiler now also exposes **[exact first-order operator tangents](reduced/operator_tangent.py)**.
Given (dG/dtheta) and (dC/dtheta), it differentiates
`G,C -> P,X -> implicit NMDA -> soma trace` analytically, reusing the local
Newton Jacobian rather than finite-differencing the solver. The compiler
matrices and complete nonlinear soma-trace tangent are checked against centered
finite differences, including a multi-parameter path. This is the first step
from geometry as fixed input to geometry as a controllable state variable; it
is not yet a claim about an optimal biological shape.

The one retained outlier is rat L6 IPC at 11.2% cell-median soma error, matching
its already-known passive graph discretization weakness.

See
**[CAUSAL_NONLINEAR_GRAPH_DIAGNOSIS.md](results/CAUSAL_NONLINEAR_GRAPH_DIAGNOSIS.md)**
and the failed precursor
**[CROSS_CELL_NONLINEAR_GRAPH_AUDIT.md](results/CROSS_CELL_NONLINEAR_GRAPH_AUDIT.md)**.

## The whole morphology graph is the cross-cell coordinate

The failed cross-cell feature regressions turned out to be asking the wrong
question.

A direct matched-passive cable solver was built from each released morphology's
**full loaded section graph** — lengths, diameters, membrane areas, topology and
the fixed passive constants — with no cross-cell fitting.

Across all **24 FCI morphologies × 6 apical branches = 144 operator packs**:

~~~text
direct morphology graph

median joint G/T NRMSE          0.21%
median local G NRMSE            0.13%
median soma T NRMSE             0.24%

median cell-level error         0.18%
cells <= 10%                    23 / 24
~~~

Classification:
`MORPHOLOGY_GRAPH_GENERATES_PASSIVE_OPERATOR`.

The extreme human L5 morphology `2057`, which the scalar morphology map missed
by nearly 8× error, is reconstructed by the direct graph at roughly **0.01%**
cell-level error.

So the cross-cell result is not "find a better morphology feature vector." It is:

~~~text
morphology graph
      |
      v
passive cable construction
      |
      +--> local Green matrix G
      +--> soma transport T
~~~

The operator becomes low-dimensional **after** the cable physics acts; forcing
the input morphology through a tiny scalar chart first destroys important
boundary/load information.

See
**[DIRECT_CABLE_GRAPH_AUDIT.md](results/DIRECT_CABLE_GRAPH_AUDIT.md)**.

## Cross-cell generalization boundary

The hard next test was **not** another within-cell perturbation. We measured
matched-passive operators on all **24 released FCI morphologies** and asked a
leave-one-cell-out model to predict an unseen neuron's local Green matrix and
site-to-soma transport from morphology alone.

~~~text
24 cells x 6 deterministic apical sections = 144 operators

training-basis PCA oracle       3.07% joint NRMSE
morphology predictor           35.22%
nearest training branch        31.99%
training mean                  66.86%

morphology beats nearest       11 / 24 cells
~~~

Classification:
`CROSS_CELL_OPERATOR_LOW_DIMENSIONAL_BUT_MORPHOLOGY_MAP_WEAK`.

So the operator family itself remains highly compressible, but the obvious
length/diameter/path/tree descriptors do **not** provide a sufficiently good
cross-cell coordinate system. Within-cell geometry interpolation survives;
arbitrary morphology-to-operator prediction does not.

See
**[CROSS_CELL_OPERATOR_AUDIT.md](results/CROSS_CELL_OPERATOR_AUDIT.md)**.

## The reduced circuit generalizes in time

Without any new fit, the same small circuit was tested on synchronous,
forward, reverse and widely staggered three-site event programs across all six
branches and all three geometries.

~~~text
72 full-model comparisons

transport oracle              0.37% median NRMSE
reduced nonlinear circuit     0.50%
reduced current waveform      0.26%

open-loop no-feedback model  41.41%
~~~

Every fixed-point solve converged and the reduced nonlinear circuit beat the
open-loop model in **72 / 72** cases.

Classification:
`TEMPORAL_GREEN_CIRCUIT_GENERALIZES_WITHOUT_REFIT`.

See
**[TEMPORAL_GREEN_CIRCUIT_AUDIT.md](results/TEMPORAL_GREEN_CIRCUIT_AUDIT.md)**.

The reduced equations and current validity envelope are written explicitly in **[MODEL.md](MODEL.md)**.

## Portable runtime

The earned reduced operator now has a **NEURON-free NumPy implementation** in
**[reduced/green_circuit_numpy.py](reduced/green_circuit_numpy.py)**.

It accepts arbitrary compatible local Green kernels, output transport kernels,
baseline trajectories and AMPA/raw-NMDA conductance programs.

**Important boundary:** that portable module currently uses the original global
waveform fixed-point algorithm. It is validated on the cell-1125 regimes above,
but the 24-cell cross-cell audit showed that this solver is not universally
robust. The cross-cell passing result uses a causal state-space graph solver
with a three-dimensional implicit solve at each time step.

~~~text
5 focused unit tests
5 passed
NEURON dependency: none
~~~

Cell-1125-specific derived kernels are intentionally not bundled with this
portable runtime while model-data provenance/licensing is kept separate from
the generic operator code. See **[reduced/README.md](reduced/README.md)**.

## The full branch now reduces to a tiny nonlinear circuit

The strongest current result removes the cached nonlinear current waveform
entirely.

For each compact three-site branch, the reduced model contains only:

~~~text
released AMPA/NMDA law
        x
3x3 local Green matrix
        x
site-to-soma transport kernels
~~~

The local currents are solved self-consistently from voltage-dependent NMDA
feedback. Across **54 comparisons** spanning six real cell-1125 branches,
original/0.8x/1.2x intrinsic length, and three input patterns:

~~~text
transport oracle soma error        0.40% median NRMSE
fully reduced circuit              0.43%
reduced current-waveform error     0.38%

held-out frozen-current model      2.82%
held-out reduced circuit           0.43%
reduced wins                       32 / 36
~~~

Classification:
`LOCAL_GREEN_MATRIX_X_SYNAPSE_LAW_REDUCES_RELEASED_NEURON`.

This is the current reusable mathematical object of Operaattori. See
**[ARCHITECTURE.md](ARCHITECTURE.md)** and
**[GREEN_CIRCUIT_AUDIT.md](results/GREEN_CIRCUIT_AUDIT.md)**.

## The architecture now composes

The latest held-out test separates the released model into two reusable pieces:

~~~text
input
  |
  v
local nonlinear current operator N_b
  |
  v
geometry-dependent transport operator T_g
  |
  v
soma
~~~

N_b was measured only in the original geometry. Branch length was then changed
to 0.80x or 1.20x and only T_g was replaced. Across 12 held-out perturbations:

~~~text
frozen original-soma attacker       9.82% median trace NRMSE
factorized T_g[N_original]          5.12%
transport oracle T_g[N_holdout]     0.46%

factorized beats frozen             9 / 12
~~~

Classification:
`TRANSPORT_X_LOCAL_NONLINEAR_OPERATOR_FACTORIZATION`.

A second attacker then reused each held-out geometry's **same transport
operator across three distinct input patterns** (middle single site, outer
pair, and triple). Across 36 holdouts:

~~~text
frozen-soma attacker       9.31% median NRMSE
factorized                 2.82%
transport oracle           0.39%
factorized wins            28 / 36
~~~

Pattern medians were 1.25%, 3.73% and 5.12% respectively.

Classification:
`TRANSPORT_OPERATOR_REUSES_ACROSS_INPUT_PATTERNS`.

This is the first result in the repo that behaves like a reusable architecture
rather than merely a causal boundary. See **[ARCHITECTURE.md](ARCHITECTURE.md)**,
**[the factorization receipt](results/OPERATOR_FACTORIZATION.md)** and
**[the cross-input attacker](results/CROSS_INPUT_TRANSPORT_AUDIT.md)**.

The repo started as a synthetic moving-operator project. From Gate 10 onward it
became a different and much more concrete investigation:

> **If a real dendritic arbor is represented as a branching field of local
> matrices, which kinds of geometric change actually change what the neuron
> computes?**

That is now the front-door question.

## The result so far: geometry is not one thing

All four rows below were measured on, or grounded in, the same reconstructed
human L2/3 pyramidal neuron, cell 1125.

| channel | intervention | measured consequence | status |
|---|---|---:|---|
| **pure 3-D embedding** | rotate one local SE(3) frame; move distal arbor **184.107 um** while preserving the cable metric | passive-transfer change **5.999e-14** | **null** |
| **intrinsic cable geometry** | stretch the same subtree metric by 20% | median passive-transfer change **36.34%** | **strong** |
| **world-space sampling, passive readout** | isometric bend changes where 1,809 material nodes sample a fixed smooth field; total drive matched | median soma change **2.11%**; 19.9% of conditions >5% | **weak / heterogeneous** |
| **nonlinear branch compartments** | coactivate compact released-model NMDA synapses inside vs across branches | median within-branch nonadditivity **70.19%**; median extra cross-branch interaction **2.01%** | **strong local modularity** |

This table is the current result.

See **[CHANNELS.md](CHANNELS.md)** for the exact measurements, commands and
stopping rules.

### The important boundary

A visibly dramatic matrix motion is not automatically a new neuronal operator.

For the ordinary cable equation,

~~~text
same length + same radius + same topology
              |
              |  arbitrary rigid re-embedding in XYZ
              v
       same intrinsic cable operator
~~~

Gate 22 moved a large real subtree almost two hundred microns and the passive
electrical model was invariant to numerical precision.

So the original "moving matrix neuron" idea has become more precise:

~~~text
local matrix motion
      |
      +--> changes intrinsic metric/radius/topology --> cable operator changes
      |
      +--> changes position in world space ----------> what the arbor encounters changes
      |
      +--> pure isometry with neither ----------------> drawing changes, cable does not
~~~

That is not new cable theory. Its value here is disciplinary: it prevents the
repo from calling every beautiful geometric deformation "computation."

The mathematical version of this boundary is in **[SYMMETRY.md](SYMMETRY.md)**.

A direct released-model audit then connected this strong metric channel to the
six Gate-20 nonlinear branches. Stretching each selected branch section by 20%
changed passive site properties by **8.77% median**, but changed the exact
single-site-subtracted nonlinear interaction by only **2.74% median**; just
**1 / 6** branches exceeded 5%. In this assay, **metric changes transport more
than the normalized nonlinear law**. See
**[METRIC_NONLINEARITY_AUDIT.md](results/METRIC_NONLINEARITY_AUDIT.md)**.

## Where the nonlinear result enters

Gates 16–20 then found something the passive picture does not contain.

Using the pinned released cell-1125 NEURON model and its human AMPA/NMDA
kinetics, compact branch coactivation survives increasingly hard attackers:
same-site single-event superposition, the paper's Hybrid-B gamma control, and a
sitewise rest-matched gamma=0.062 control.

Gate 20's hierarchical test is the cleanest summary:

~~~text
HUMAN
median within-branch nonlinearity        70.19%
median extra cross-branch interaction     2.01%
median modularity margin log              0.4582
positive modularity pairs                14 / 15

rest-matched gamma=0.062 margin           0.0579
HUMAN extra gamma-specific margin         0.3971
~~~

Classification:

~~~text
HUMAN_GAMMA_STRENGTHENS_SEMI_INDEPENDENT_COMPARTMENTS
~~~

That earns a useful abstraction:

> **one real dendritic scaffold contains multiple semi-independent nonlinear
> subunits.**

Gate 21 then showed that redistributing the same 48-synapse budget across those
subunits changes somatic output strongly, and the nonlinear branch basis
predicts which redistribution wins much better than independent sites. But it
does **not** reconstruct the complete soma waveform closely enough to call the
whole neuron a clean sum of branch transfer functions.

## The matrix-exponential connection

The mathematical seed for this line is:

~~~text
dx/dt = A x
x(t)  = exp(t A) x(0)
~~~

A matrix generator becomes continuous motion: rotation, stretch, shear, spiral.

Operaattori does **not** simply insert an arbitrary 3x3 flow into a neuron and
call the result biology. The real morphology is encoded as parent-local 4x4
transforms,

~~~text
W_child = W_parent T_child
~~~

so a local transform can move an entire distal subtree coherently. Gate 10
shows the representation; Gate 22 tells us which components of that motion are
physically causal for an intrinsic cable model.

The old synthetic Lie-generator/moving-operator work is preserved in
**[PAST.md](PAST.md)**. It is provenance, not the current headline.

## Gate 23 is a toy spatial-sampling assay, not extracellular-field physics

Gate 23 deliberately used a normalized positive exponential world-space drive.
It showed that relative scaffold/world geometry can be causally separated from
the intrinsic cable operator, but its 2.11% median effect is weak.

It should **not** be reclassified as a reproduction of the known bent-cable
projection formula. That literature uses extracellular potential in the cable
equation.

Relevant prior work:

- Aspart, Remme & Obermayer (2018), *Differential polarization of cortical
  pyramidal neuron dendrites through weak extracellular fields*:
  https://doi.org/10.1371/journal.pcbi.1006124
- Fan et al. (2023/2024), *Electric field effects on neuronal input-output
  relationship by regulating NMDA spikes*:
  https://doi.org/10.1007/s11571-022-09922-y

The first paper derives the bent-cable projection geometry and reports
frequency-dependent dendritic polarization. The second explicitly shows that
electric-field polarization can shift NMDA-spike generation.

Therefore Operaattori must **not** claim that "field + NMDA" is an undiscovered
phenomenon.

## Gate 24 closes the world-field loop — but the nonlinear effect is tiny

Gate 24 replaced Gate 23's toy current redistribution with a real uniform
extracellular potential coupled through NEURON while keeping the intrinsic
cell-1125 cable unchanged.

The controls are exact:

~~~text
zero-field embedding difference          0
material-locked difference               0
~~~

Under the locked 1 V/m, 35-degree, 24-synapse protocol:

~~~text
proximal HUMAN bend factor          1.0003x
distal HUMAN bend factor            1.0002x
frequencies with >5% proximal effect   0 / 3
~~~

Classification:

~~~text
FIELD_COUPLING_PRESENT_BUT_NOT_COMPARTMENT_SELECTIVE
~~~

So a real world-space field gives the isometric embedding a physical coupling,
but under this weak-field protocol it does **not** materially change the
already-measured nonlinear branch interaction.

A bounded follow-up amplitude audit then held branch, bend, synaptic dose and
frequency fixed and varied only the field from 0.25 to 2 V/m. The bend effect
was almost perfectly linear through the origin (**R² = 0.999970**), with a
first-order scale of about **175 V/m** for the old 5% effect ruler. That value
is an extrapolated scale, not a simulated high-field result. See
**[FIELD_SCALING_AUDIT.md](results/FIELD_SCALING_AUDIT.md)**.

This closes rather than rescues the weak-field route.

No Gate 25 is open.

## Scaffold-line gate map

| gate | question | result |
|---|---|---|
| [10](results/GATE10.md) | can cell 1125 be represented exactly as parent-local matrices? | yes |
| [11](results/GATE11.md) | can heterogeneous cable matrices be order-sensitive? | yes |
| [12](results/GATE12.md) | is the biological fine serial order special beyond taper? | no; monotonic taper explains it |
| [13](results/GATE13.md) | do branch junctions matter beyond the same serial path? | yes; strong branch loading |
| [14](results/GATE14.md) | does reduced human NMDA add branch nonlinearity? | modest gain, no new dimension |
| [15](results/GATE15.md) | does quasi-static compact locality selectively amplify human NMDA? | no |
| [16](results/GATE16.md) | does the released time-domain model change that verdict? | yes |
| [17](results/GATE17.md) | does it survive exact same-site superposition? | yes, heterogeneous |
| [18](results/GATE18.md) | does it survive Hybrid-B / rest-matched gamma? | narrow yes |
| [19](results/GATE19.md) | is a compact branch a 4-ms sequence-direction processor? | no |
| [20](results/GATE20.md) | are nonlinear interactions stronger within than between branches? | **yes** |
| [21](results/GATE21.md) | do branch subunits cleanly factorize fixed-budget soma traces? | mixed; scalar signature yes, full trace no |
| [22](results/GATE22.md) | does pure 3-D re-embedding change cable computation? | **no, to numerical precision** |
| [23](results/GATE23.md) | does passive toy world sampling make that embedding robustly functional? | weak |
| [24](results/GATE24.md) | real extracellular field x movable scaffold x nonlinear compartments | field coupling present; **nonlinear bend effect ~0.03%, not selective** |

## Reproduce the current headline rows

Gate 20 and Gate 21 require the pinned released FCI model and NEURON runtime;
the CI workflow shows the complete environment setup.

~~~bash
# Pure embedding vs intrinsic metric
python experiments/gate22_embedding_vs_metric.py \
  --paths 32 --twist-deg 35 --stretch 1.20

# Passive toy world-space sampling
python experiments/gate23_spatial_field.py \
  --directions 48 --twist-deg 35
~~~

Released-model commands are documented in [CHANNELS.md](CHANNELS.md) and
`.github/workflows/ci.yml`.

## What this repo is not claiming

It is not claiming:

- that matrix exponentials are new;
- that passive cable invariance under rigid re-embedding is new physics;
- that extracellular fields modulating NMDA spikes are new;
- that cell 1125 is "intelligent";
- that the scaffold already grows usefully.

What the repo **does** have is a progressively attacked causal decomposition of
a moving real-neuron scaffold, with positive and negative results kept in the
ledger.

That is the thing to preserve.
