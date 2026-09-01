# Metric → nonlinear branch audit

The weak-field embedding route is now closed at the audited scale. The strong
unconnected arrow is intrinsic metric → nonlinear branch computation.

Gate 22 showed that changing intrinsic cable length strongly changes passive
transfer. Gate 20 independently established six compact nonlinear branch
subunits. This audit asks whether those two results actually meet.

## Locked intervention

Use the same pinned FCI human cell-1125 model and the exact Gate-20 six-branch
basis.

For each of the six branches:

1. measure the HUMAN three-site nonlinear interaction in the untouched model;
2. build a fresh identical model;
3. increase **only that selected section's cable length by 20%**;
4. keep diameter, topology, section identity, normalized synapse locations,
   HUMAN NMDA kinetics and multiplicity unchanged;
5. recompute the exact same interaction and passive site features.

The active dose is locked to Gate 20:

~~~text
3 sites x multiplicity 8 = 24 simultaneous virtual synapses
~~~

No field and no XYZ re-embedding are involved.

This is a causal geometry intervention, not a claim that a living dendrite
acutely stretches by 20%.

## Nonlinear ruler

For each geometry,

~~~text
I = AUC(simultaneous three-site local depolarization)
    / AUC(sum of the three exact single-site depolarizations)
~~~

The branch-level metric effect is

~~~text
B = abs(log(I_stretched / I_original))
~~~

The audit earns

~~~text
INTRINSIC_METRIC_MODULATES_NONLINEAR_BRANCH_INTERACTION
~~~

only if:

- median B >= log(1.05);
- at least 4/6 branches individually exceed log(1.05);
- the passive positive control changes site input/transfer by >=1% median;
- no simultaneous HUMAN run reaches the conservative -20 mV soma guard.

Otherwise the result is

~~~text
METRIC_CHANGES_PASSIVE_TRANSPORT_BUT_NOT_NONLINEAR_RATIO
~~~

provided the passive positive control passes.

The point is not to demand a positive result. The interaction ratio already
subtracts the exact single-site responses, so failure would mean the metric
edit changes transport without substantially changing the branch's nonlinear
superposition law.

This is a focused audit, not Gate 25.


## Receipt

The corrected locked audit completed successfully. The first CI attempt had
only a template-instance naming bug (`Model[0].apic[100]` versus
`Model[1].apic[100]`) and produced no scientific result. The rerun compared
the stable section identities without changing any experimental parameter.

~~~text
branch       length edit (um)       passive factor   I original -> stretched   nonlinear factor
apic[100]    385.4 -> 462.5         1.0928x          1.2688 -> 1.2474           1.0172x
apic[77]     359.2 -> 431.1         1.0845x          1.2383 -> 1.3120           1.0595x
apic[96]     353.9 -> 424.7         1.1014x          2.0496 -> 2.0230           1.0131x
apic[64]     324.2 -> 389.0         1.0815x          2.0796 -> 2.1198           1.0193x
apic[58]     316.9 -> 380.2         1.0909x          1.9727 -> 2.0639           1.0462x
apic[69]     316.5 -> 379.7         1.0770x          1.4683 -> 1.5205           1.0355x

aggregate
  median passive site factor                      1.0877x
  passive branches >1%                            6 / 6
  median nonlinear interaction factor             1.0274x
  nonlinear branches >5%                          1 / 6
  spike guard fraction                            0
~~~

Classification:

~~~text
METRIC_CHANGES_PASSIVE_TRANSPORT_BUT_NOT_NONLINEAR_RATIO
~~~

The intrinsic intervention is electrically real in every branch, but the exact
three-site interaction ratio is much more stable than the passive site
properties. A 20% local length edit therefore does not robustly turn into a
change in the branch's nonlinear superposition law under the locked Gate-20
dose.

This suggests a useful separation:

> **intrinsic metric strongly controls transport, while the local nonlinear
> interaction law can remain comparatively stable after exact single-site
> normalization.**

That is not the same as saying morphology is irrelevant to NMDA responses.
Absolute single-site and somatic responses can change with cable geometry, and
prior literature already establishes morphology-dependent NMDA-spike
thresholds. The narrower Operaattori result concerns the normalized nonlinear
interaction law on these six already-established cell-1125 compact branches.

No stretch or dose rescue scan is opened.
