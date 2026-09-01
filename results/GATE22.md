# Gate 22 — which geometry in the matrix scaffold is electrically causal?

Gate 21 did not earn a clean branch-transfer factorization, so this gate does
**not** rescue it with new doses. Instead it returns to the original scaffold
question.

Gate 10 showed that one local SE(3) rotation can visibly bend a large distal
subtree while preserving every cable length. The important question is now:

> Does that change the intrinsic cable operator merely because the neuron looks
> different in 3-D?

For the ordinary passive cable equation the answer should be **no**. Internal
transport depends on the neurite metric and radius, not on how an isometric
copy of the same cable is embedded in external XYZ space.

That is a crucial boundary for the geometric-neuron idea.

## Locked comparison

Use the pinned cell-1125 reconstruction and Gate-10 matrix scaffold.

Choose the same deterministic internal bifurcation rule and compare:

1. **ORIGINAL** — untouched scaffold.
2. **SE3 TWIST** — rotate one parent-local frame by 35 degrees. The complete
   distal subtree moves in 3-D, but local translations, cable lengths, radii
   and topology are preserved.
3. **METRIC STRETCH** — leave the attachment and topology fixed but multiply
   every parent-local translation wholly inside that same subtree by 1.20.
   Radii are unchanged. This changes cable length and therefore the intrinsic
   electrical metric.

No learning and no synaptic nonlinearity are involved.

## Electrical ruler

For up to 32 affected dendritic root-to-tip paths, evaluate the same classical
passive ABCD cable matrices used in Gate 11 over 1–300 Hz.

Report changes in:

- sealed input impedance;
- sealed distal voltage gain.

## Preregistered interpretation

The intended boundary is earned only if all of these hold:

- the SE3 twist moves the distal arbor by more than 10 um;
- its maximum cable-length change is below 1e-7 um;
- its maximum passive-transfer difference is below 1e-9;
- the 20% metric stretch causes at least a 1% median transfer change.

If so the classification is:

~~~text
CABLE_MODEL_IGNORES_ISOMETRIC_3D_EMBEDDING
~~~

This does **not** say real neurons ignore 3-D position. Spatial embedding can
matter through extracellular fields, spatially located synapses, diffusion,
contact, mechanics and other external couplings. It says something narrower
and important:

> an SE(3) scaffold rotation is not automatically a new internal electrical
> operator when length, radius and topology are unchanged.

## Stopping line

If this boundary passes, the next honest route for 3-D scaffold geometry is an
**extrinsic spatial coupling** test: put the same neuron in a fixed structured
external field and ask whether an isometric local scaffold bend changes what
the branches sample. Do not pretend pure XYZ bending changes cable physics by
itself.


## Receipt

The locked assay passed cleanly on the real cell-1125 scaffold:

~~~text
affected dendritic paths                13
pivot node                            9418

35-degree SE3 bend
  max distal 3-D displacement       184.107 um
  max cable-length change          1.483e-13 um
  max passive-transfer change      5.999e-14

20% intrinsic metric stretch
  median passive-transfer change      0.3634
  paths with >1% change              13 / 13
~~~

Classification:

~~~text
CABLE_MODEL_IGNORES_ISOMETRIC_3D_EMBEDDING
~~~

The same local matrix operation can therefore be visually enormous and
electrically null. In the classical internal cable model, what matters is the
intrinsic neurite metric/radius/topology. Absolute XYZ embedding is not an
independent electrical degree of freedom.

This is a boundary on the scaffold idea, not a claim that biological neurons
ignore space. It tells us exactly what must be added for a pure bend to matter:
some process that lives in world space rather than only along the cable.
