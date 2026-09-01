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
