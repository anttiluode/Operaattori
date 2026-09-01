# Operaattori

![cell 1125 reconstructed as a local-matrix scaffold](cell1125_original_scaffold.png)

**Operaattori** is Finnish for **operator**.

## [▶ Live Bend Fence demo](https://anttiluode.github.io/Operaattori/)

**Drag a real human dendritic subtree through space.** A 35° rigid bend moves it
**184.107 µm** while the classical passive cable transfer changes by only
**5.999e-14**. Switch to intrinsic stretch and a much less dramatic-looking
20% metric change produces a **36.34% median passive-transfer change**.

The browser computes the geometry live and displays the locked Gate-22
electrical measurements as measured anchors; it does not fabricate interpolated
electrophysiology.

A separate **72-bend symmetry audit** then repeated the null across six
bifurcations, all three local axes and positive/negative bend angles: the
largest subtree motion was **619.350 µm**, while the largest passive-transfer
change was only **9.710e-14**. The 20% intrinsic-metric controls changed
transfer by **11.84% median**. See
**[results/SYMMETRY_AUDIT.md](results/SYMMETRY_AUDIT.md)**.

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
