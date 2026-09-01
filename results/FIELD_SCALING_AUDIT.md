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


## Receipt

The bounded audit passed in its local-linearity sense.

~~~text
proximal section                      apic[77]
max bend displacement                 278.569 um
frequency                             15 Hz

field        signed log bend effect       factor
0.25 V/m     -6.856178e-05                1.000069x
0.50 V/m     -1.375044e-04                1.000138x
1.00 V/m     -2.765310e-04                1.000277x
2.00 V/m     -5.591471e-04                1.000559x

through-origin slope                  -2.787234e-04 log/(V/m)
through-origin R2                      0.999970
sign consistent                       yes
spike guard fraction                  0
local first-order field for 5%        175.049 V/m
~~~

Classification:

~~~text
WEAK_FIELD_BEND_EFFECT_IN_LOCAL_LINEAR_REGIME
~~~

Within 0.25–2 V/m, the geometry-dependent nonlinear interaction perturbation is
therefore an extremely clean small-signal response. Doubling field amplitude
approximately doubles the signed effect.

The 175.049 V/m value is **not** a simulation result at 175 V/m and is not a
claim that the response stays linear there. It is the first-order scale implied
by the measured weak-field slope.

For context, conventional human transcranial electric stimulation generally
produces intracranial fields below about 1 V/m. Huang et al. (2017) measured and
modeled cortical maxima around 0.4 V/m per 1 mA and about 0.8 V/m for typical
2 mA stimulation:

https://doi.org/10.7554/eLife.18834

A review of immediate TES physiology likewise summarizes conventional human
fields as <1 V/m:

https://doi.org/10.1038/s41467-018-07233-7

So the honest conclusion is not "try 175 V/m." It is:

> **The weak-field isometric-bend route is physically real but far too small,
> in this locked cell-1125 assay, to be the nonlinear computational amplifier
> we were looking for.**

No larger-amplitude rescue scan is opened.
