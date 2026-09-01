# Embedding symmetry of the scaffold

Operaattori began with moving matrices. Gates 22–24 make the physically useful
distinction sharper:

> **A local matrix motion matters only through the physical quantities it
> changes.**

For the classical one-dimensional cable model, an important class of matrix
motions changes no intrinsic electrical quantity at all.

## 1. Intrinsic cable operator

Write one dendritic branch using material arc length (s). In schematic form a
passive cable equation is

[
c_m(s),partial_t v
=
partial_s!left(g_a(s),partial_s vight)
-
g_m(s),[v-E_L]
+
i_{mathrm{material}}(s,t).
]

The exact coefficient conventions are not important here. What matters is what
the operator sees:

- material arc length;
- radius / cross-section;
- membrane and axial parameters;
- branch connectivity;
- material-locked inputs and boundary conditions.

The absolute world-space coordinate

[
mathbf r(s)inmathbb R^3
]

does not appear in the ordinary intrinsic cable operator.

Therefore, if a new embedding (mathbf r'(s)) preserves the material metric,
radius, topology and material inputs, the intrinsic operator is unchanged even
if the dendrite looks dramatically different in XYZ.

That is the symmetry tested visually by Gate 22 and the Bend Fence demo.

## 2. This invariance is known cable geometry, not an Operaattori novelty claim

López-Sánchez & Romero derived a generalized cable equation for curved cable
geometry using a Frenet–Serret description:

> *Cable equation for general geometry*, Physical Review E 95, 022403 (2017)

https://doi.org/10.1103/PhysRevE.95.022403

For a constant circular cross-section, their generalized equation depends on
neither curvature nor torsion.

Operaattori's contribution is therefore not "curvature does not matter." The
repo turns that known boundary into an explicit causal control on a real human
reconstruction represented as parent-local matrices.

## 3. What breaks the symmetry?

A world-space field introduces a function of position, for example an
extracellular potential

[
V_e(mathbf r,t).
]

The cable now samples the pullback

[
V_e(mathbf r(s),t).
]

After an isometric re-embedding,

[
mathbf r(s)ightarrowmathbf r'(s),
]

the intrinsic cable is unchanged, but generally

[
V_e(mathbf r(s),t)
eq V_e(mathbf r'(s),t).
]

So the *relationship between scaffold and world* becomes causal.

The same pattern applies to other possible external couplings:

- spatially located synaptic partners;
- extracellular electric fields;
- diffusion or chemical gradients;
- collisions/contact;
- mechanical forces;
- tissue boundaries.

These are not automatically important. They merely provide a channel through
which embedding can matter.

Gate 23 added an abstract spatial-sampling channel and found a weak,
heterogeneous effect. Gate 24 used a real extracellular potential and found
physical coupling but only about a 0.03% change in the locked nonlinear branch
interaction.

## 4. What breaks the intrinsic symmetry directly?

Change the metric itself.

Gate 22 multiplies parent-local translations wholly inside one subtree by 1.20.
That changes cable lengths while preserving radius and topology.

The result:

~~~text
rigid local bend
  maximum XYZ displacement          184.107 um
  maximum passive-transfer change   5.999e-14

20% intrinsic metric stretch
  median passive-transfer change    36.34%
~~~

The visually larger operation is electrically null; the visually smaller
metric operation is electrically strong.

## 5. Matrix interpretation

The morphology is reconstructed from parent-local homogeneous transforms,

[
W_i = W_{p(i)}T_i.
]

A change to the rotational block of one (T_i) can re-embed a complete distal
subtree while leaving its material translations unchanged.

A change to the translation magnitudes changes the intrinsic metric.

So the same matrix scaffold carries physically different channels:

~~~text
local matrix edit
      |
      +-- rotation only, metric preserved
      |       |
      |       +-- no world coupling --> intrinsic cable symmetry
      |       |
      |       +-- world coupling ----> embedding can become causal
      |
      +-- translation magnitude / radius / topology changed
              |
              +----------------------> intrinsic cable operator changes
~~~

This is the central fence of the current repo.

## 6. Why the symmetry audit exists

Gate 22 used one selected bifurcation and one 35-degree bend. Because the live
demo now makes that example prominent, the implementation should not rely on
one lucky geometry.

`experiments/symmetry_audit.py` tests multiple bifurcations, all three local
axes and positive/negative bend angles. The audit is expected to reproduce the
known invariance. Its purpose is implementation integrity, not discovery.

See [results/SYMMETRY_AUDIT.md](results/SYMMETRY_AUDIT.md).
