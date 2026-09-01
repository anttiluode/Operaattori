# Gate 23 — can a fixed world field make an isometric scaffold bend functional?

Gate 22 established a hard boundary:

    same cable metric + same radii + same topology
    + very different XYZ embedding
    -> same intrinsic passive cable operator

So Gate 23 adds the smallest missing ingredient: a **fixed spatial environment**.

This is deliberately an abstract spatial-drive assay, not a claim that the
chosen field is a literal extracellular or synaptic mechanism.

## Object

Use the same cell-1125 matrix scaffold and the same deterministic Gate-10/22
pivot.

Compare the original arbor and the 35-degree isometric bend from Gate 22.

The intrinsic electrical tree is kept exactly the same. Only the material
nodes' world coordinates change.

## Fixed world fields

Across the moved subtree, define smooth positive spatial gradients

    phi_d,beta(x) = exp(beta * d dot (x - x_pivot) / R)

for 48 deterministic directions d on the sphere and three locked gradient
strengths:

    beta = 0.5, 1.0, 2.0

R is the RMS spatial radius of the original moved subtree around the pivot.

Each cable node samples the field at its own world coordinate. The sampled
drive is weighted by local cable surface proxy 2*pi*r*L.

### Equal-total-input attacker

For every field separately, normalize ORIGINAL and BENT injections to the same
total current.

Therefore a scalar ruler that knows only total exposure predicts **no output
change**.

Any remaining change comes from redistributing the same total drive among
material nodes with different electrical transfer to the soma.

## Electrical readout

Use the exact whole-tree passive transfer coefficients from Gate 13 at:

    5, 20, 80 Hz

For an injection pattern I_i,

    y_soma = sum_i T_i I_i

where T_i is unchanged by the isometric bend.

## Material-locked control

As a sanity attacker, keep each material node's ORIGINAL sampled field value
attached to that node even after the geometric bend.

Because neither T_i nor I_i then changes, the somatic response must be
identical to numerical precision.

This separates:

    changing the intrinsic operator

from

    changing what the same operator samples in world space.

## Preregistered positive classification

The result earns

~~~text
SPATIAL_COUPLING_MAKES_EMBEDDING_FUNCTIONAL
~~~

only if:

- the isometric bend still changes intrinsic passive transfer by <1e-9
  (Gate-22 invariant);
- material-locked output difference is <1e-12;
- fixed-world, equal-total-input median output difference is >=5%;
- at least half of field/frequency conditions exceed 5%.

A pass earns only this statement:

> Isometric scaffold geometry can become functionally causal when material
> locations couple to a structured world-space field.

It does not earn growth, learning, biological optimality or intelligence.

## Stopping line

If Gate 23 passes, the first defensible growth question is no longer “can a
branch bend?” It is:

> can a local developmental rule change intrinsic or extrinsic geometry so
> that later sampling/performance improves under a fixed environment?

That is the point where adaptive scaffold motion becomes worth testing.


## Receipt — causal coupling exists, but the locked smooth-field effect is weak

The 35-degree isometric bend was evaluated across all 48 field directions,
three gradient strengths and three frequencies:

~~~text
sampled material nodes                    1809
field/frequency conditions                 432

max intrinsic cable-length change      1.426e-13 um
material-locked max output change       0.000e+00

fixed-world equal-total condition
  median output change                     0.0211
  conditions >5%                            0.199
  conditions >10%                           0.060
  maximum output change                     0.2230
~~~

Classification:

~~~text
FIXED_WORLD_FIELD_EFFECT_WEAK
~~~

The zero material-locked error verifies the causal decomposition: the bend does
not secretly alter the intrinsic passive operator. Fixed-world fields can
produce sizeable effects for some orientations, but the preregistered family is
not robust enough to earn the stronger claim that spatial coupling generally
makes this passive embedding functional.

We do not rescue the gate by sharpening the field, changing the normalization
or selecting favorable directions after seeing the result.

The combined boundary after Gates 20–23 is now informative:

- real branches can be strong local nonlinear compartments;
- a pure SE(3) bend is intrinsically cable-invariant;
- smooth spatial resampling through a passive soma is mostly weak.

A future Gate 24 should therefore not be “more gradient.” The clean candidate
is to combine two mechanisms that independently survived earlier attacks:
**world-space sampling and the measured nonlinear branch compartments**. That
would require a fresh locked protocol and a real nonlinear model, not a tuned
version of Gate 23.
