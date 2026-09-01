# Symmetry audit — is Gate 22 one lucky bend?

This is **not Gate 25** and it does not search for a larger effect.

Gate 22 used one deterministic 35-degree local bend on cell 1125 and found a
visually large but electrically null intervention in the classical passive
cable model. The live Bend Fence demo now makes that result prominent, so the
implementation owes a broader audit:

> does the same intrinsic-cable invariance survive many local SE(3)
> re-embeddings at different bifurcations, axes and signs?

## Literature fence

This invariance is not being claimed as new theory.

López-Sánchez & Romero, *Cable equation for general geometry*, Physical Review
E 95, 022403 (2017), derive a generalized cable equation using a Frenet-Serret
description. For constant circular cross-section, their equation depends on
neither curvature nor torsion.

https://doi.org/10.1103/PhysRevE.95.022403

So the purpose here is narrower: verify that Operaattori's real-neuron matrix
scaffold and passive transfer implementation actually respect that known
symmetry rather than producing Gate 22 because of one special pivot.

## Locked audit

Use the same pinned human L2/3 cell-1125 reconstruction.

Select six deterministic non-axonal bifurcations spanning viable subtree sizes.
At each pivot test:

- local axes x, y and z;
- angles -35, -15, +15 and +35 degrees;
- three affected dendritic tip paths;
- passive transfer at 1, 5, 20, 80 and 300 Hz.

That gives 72 independent local rigid-bend trials.

For every bend, measure:

- maximum distal XYZ displacement;
- maximum cable-length change;
- maximum passive-transfer change.

At each of the same pivots, run a 20% intrinsic subtree stretch as the positive
control.

## Locked rulers

The audit passes only if:

~~~text
max rigid-bend cable-length change          < 1e-7 um
max rigid-bend passive-transfer change      < 1e-9
at least one visual displacement            > 10 um
median metric-control transfer change       >= 1%
fraction metric controls >1%                >= 0.80
~~~

Passing classification:

~~~text
KNOWN_REEMBEDDING_SYMMETRY_REPLICATED_ON_REAL_SCAFFOLD
~~~

Failure is a software/scaffold problem to investigate. It is not permission to
retune the biology.

The result will be written to
`results/symmetry_audit/symmetry_audit.json`.


## Receipt

The locked audit passed in CI.

~~~text
bifurcations audited                    6
local rigid-bend trials                72
axes                                   x, y, z
angles                                 -35, -15, +15, +35 deg

maximum visual displacement            619.350 um
maximum rigid-bend cable-length error  1.652e-13 um
maximum passive-transfer change        9.710e-14

20% intrinsic metric controls
  median transfer change               0.1184
  controls >1%                         6 / 6
~~~

Classification:

~~~text
KNOWN_REEMBEDDING_SYMMETRY_REPLICATED_ON_REAL_SCAFFOLD
~~~

Gate 22 was therefore not a lucky pivot or axis. Across a deliberately varied
family of local SE(3) edits, the real scaffold can move hundreds of microns
while the intrinsic passive cable response remains invariant to numerical
precision. The same pivots respond immediately when their material metric is
changed.

This is a validation of the scaffold implementation against known cable
geometry, not a novelty claim about curvature or torsion.
