# Field-scaling audit — how dead is the weak-field route?

Gate 24 found real extracellular-field coupling but only a tiny change in the
nonlinear compact-branch interaction after a 35-degree isometric re-embedding.

This audit does **not** increase field strength until a desired effect appears.
It asks one bounded small-signal question:

> Is the Gate-24 bend effect locally proportional to field amplitude between
> 0.25 and 2 V/m?

If yes, the measured slope can be used for a transparent local extrapolation to
the old 5% Gate-24 ruler. If no, no such extrapolation is reported.

## Locked object

Reuse exactly the Gate-24:

- pinned FCI cell-1125 model and commit;
- proximal large-subtree arm;
- same three compact sites;
- HUMAN gamma and NMDA ratio;
- multiplicity 8 at each site = 24 virtual synapses;
- 35-degree isometric bend;
- apical principal-axis uniform extracellular field;
- 15 Hz, the center of the preregistered 10/15/20 Hz panel.

The only independent variable is field amplitude:

~~~text
0.25, 0.50, 1.00, 2.00 V/m
~~~

No branch, bend, synaptic-dose or frequency scan is allowed.

## Measurement

At each amplitude compute the same Gate-24 local nonlinear interaction ratio

~~~text
I = AUC(simultaneous local synaptic perturbation)
    / AUC(sum of exact single-site perturbations)
~~~

for ORIGINAL and BENT embeddings.

Use the signed bend effect

~~~text
y(E) = log(I_bent / I_original)
~~~

and fit a line through the physical origin:

~~~text
y(E) = slope * E
~~~

Report:

- the four signed and absolute effects;
- through-origin R²;
- sign consistency;
- spike guard;
- the local extrapolated field required to reach
  `abs(log(I_bent/I_original)) = log(1.05)`.

The extrapolated field is reported only if:

~~~text
through-origin R² >= 0.98
effect sign is consistent across all nonzero amplitudes
no run trips the somatic spike guard
~~~

## Interpretation fence

A successful linear fit does **not** establish behavior at the extrapolated
field. It says only that, within 0.25–2 V/m, the perturbation behaves like a
small-signal response and gives a first-order scale estimate.

If the estimated 5% field lies far beyond the audited panel, that closes the
weak-field route more strongly than simply saying Gate 24 was small.

This is an audit, not Gate 25.
