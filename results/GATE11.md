# Gate 11 — transport phase lives on the real matrix scaffold

Gate 10 established the geometry: the Aizenbud human L2/3 exemplar can be
represented exactly by a branching field of local SE(3) frames.

Gate 11 attaches a second matrix to the physical segments, but deliberately
uses a **classical passive-cable operator** rather than an invented neural
operator.

## Local operator

For a cylindrical segment at angular frequency `omega`:

```text
r_a = R_a / (pi a^2)
y_m = 2 pi a (1/R_m + i omega C_m)

gamma = sqrt(r_a y_m)
Z0    = sqrt(r_a / y_m)
```

and the segment's two-port matrix is

```text
          [ cosh(gamma l)       Z0 sinh(gamma l) ]
M(l,a) =  [                                         ]
          [ sinh(gamma l)/Z0    cosh(gamma l)    ]
```

mapping distal voltage/current to proximal voltage/current.

The passive parameters are the same values used for all detailed neuron models
in Aizenbud et al.:

```text
Cm = 1 uF/cm^2
Ra = 150 ohm cm
Rm = 20,000 ohm cm^2
```

This is not the paper's full neuron: Gate 11 omits NMDA, active soma/axon
channels and branch loading. It is a controlled path-level transport assay.

## Why this is the right first operator

At fixed frequency, if every segment has the same radius, all segment matrices
belong to the same one-parameter cable flow:

```text
M(l1) M(l2) = M(l1+l2) = M(l2) M(l1)
```

Order disappears.

Real dendrites change radius. Then `gamma` and `Z0` change from segment to
segment and generally

```text
M_i M_j != M_j M_i
```

So the real morphology gives Gate 9's noncommuting-path mechanism a physical
meaning: taper/diameter changes alter the transport operator.

## Assay

For the longest dendritic root-to-tip paths in cell 1125:

1. derive every segment's real length and mean radius;
2. compute its cable matrix from 1 to 300 Hz;
3. compose matrices in the real proximal->distal order;
4. reverse the exact same segment multiset;
5. randomly shuffle the exact same segment multiset;
6. create an area-matched **uniform-radius** path with the same segment lengths.

Read:

- sealed-end input impedance;
- distal voltage gain;
- gain phase;
- numerical group delay;
- local adjacent commutator action.

The reversal/shuffle controls preserve every segment. Only order changes.

The uniform-radius control preserves total path length and membrane surface
area but removes operator heterogeneity. It should be order-invariant to
roundoff.

## First full-cell attempt — functional gate positive, microscopic meter wrong

The first 64-path CI run produced a strong path-order effect but failed the
preregistered adjacent-point commutator threshold:

```text
paths audited                         64
median path length                    357.6 um
median point segments/path            241

point-adjacent commutator action      4.145e-17   <- failed meter

real vs reversed impedance            0.315524
real vs shuffled impedance            0.174423
real vs reversed gain                 0.315524
median gain phase difference          0.237898 rad
median group-delay difference         1.52477 ms

uniform-radius reverse control        1.077e-14
```

The failure is kept rather than relabeled as a pass.

The ASC morphology is sampled at roughly micron-scale point spacing. A
commutator between two immediately adjacent infinitesimal cable transfer
pieces is second-order in those tiny lengths, so that particular diagnostic
falls at numerical precision even though hundreds of pieces compose into a
large order-dependent transfer.

The gate was therefore **not weakened on the functional result**. The revised
diagnostic composes contiguous pieces into 25-um physical cable blocks and
measures block commutators. Reversal and shuffle of the exact same segment
multiset remain the primary causal tests; the area-matched uniform-radius path
remains the commuting attacker.

## What a pass means

A pass earns:

> **spatial ordering of local passive transport matrices is a real functional
> degree of freedom on the reconstructed human dendrite.**

It does **not** earn intelligence, learning, or a claim that neurons exploit
Lie brackets explicitly.

The next gate can finally ask a functional/developmental question: can a local
update law alter operator attachment or geometry so that a temporal task
improves, against shuffled-credit, commuting, fixed-scaffold and dense-matrix
attackers?
