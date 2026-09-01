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

These rows should not be merged into one vague word such as "geometry."

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
       branch compartments -------------> candidate thresholded amplification
~~~

The last arrow is not yet an Operaattori result. It is Gate 24.

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
