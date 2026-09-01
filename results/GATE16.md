# Gate 16 — last rescue: the released kinetic model

Gate 15 found strong electrical locality but no quasi-static NMDA-specific
same-branch advantage. Gate 16 tests the only remaining obvious objection:

> Did the equilibrium reduction erase a transient effect produced by the
> released AMPA/NMDA kinetics?

This gate runs the authors' pinned cell-1125 NEURON model directly. It does
not use Operaattori's passive cable solver for the dynamics.

## Matching inside NEURON

After the released model settles, NEURON's own DC Impedance calculation gives
for every dendritic segment:

- input impedance;
- transfer impedance to the soma.

For each long section, three adjacent midpoint segments fitting within a
55-um window form the compact set. They are greedily matched to three unique
segments on three different sections using log input impedance and log soma
transfer. Selection occurs before any nonlinear response is seen.

## Input

One synchronous released-model event is delivered at 60 ms. The existing
AMPANMDA_EMS NetCon weight is multiplied by 1, 4, or 16, which is exactly
equivalent to 1, 4, or 16 simultaneous identical events because NET_RECEIVE
adds the conductance state linearly.

With three sites this gives 3, 12, and 48 simultaneous virtual synapses.

## Human versus frozen block

HUMAN uses the released gamma=0.078 and human NMDA ratio.

The frozen control does not replace the mechanism. It runs the same released
AMPANMDA_EMS point process with gamma=0, making the magnesium gate constant,
and rescales NMDA_ratio so its effective NMDA conductance exactly equals the
human mechanism's block at -70 mV. Thus the control removes voltage feedback
while matching the resting effective NMDA strength.

## Primary readout

For each arrangement, integrate positive local depolarization for 90 ms after
the event. Compute:

    Rcluster = AUC(HUMAN clustered) / AUC(FROZEN clustered)
    Rspread  = AUC(HUMAN dispersed) / AUC(FROZEN dispersed)

and the temporal locality index:

    Ldynamic = Rcluster / Rspread

Soma-voltage AUC and somatic spike counts are secondary guards.

## Outcome

- DYNAMIC_NMDA_LOCALITY_ADVANTAGE_PRESENT
- NO_DYNAMIC_NMDA_LOCALITY_ADVANTAGE
- DYNAMIC_PASSIVE_MATCH_INADEQUATE
- DYNAMIC_ASSAY_SPIKING_CONFOUNDED

Gate 16 is the last allowed rescue of the Gate-15 locality hypothesis. A
negative result closes this mechanistic branch before any developmental
growth is introduced.


## Receipt — the time-domain model reverses Gate 15

The pinned FCI cell-1125 model compiled and ran under NEURON in CI. The
single-synapse smoke test produced finite AMPA and NMDA currents plus local and
somatic voltage transients using the released human parameters.

The scientific six-branch assay then gave:

~~~text
median compact span                    41.92 um
median passive Rinput match             1.037x
median passive soma-transfer match      1.126x

48-synapse local-AUC locality           1.5143
fraction branches locality >1.05        0.833  (5 / 6)
median soma-AUC locality                1.1359
branches with somatic spikes            0 / 6
~~~

Per-branch high-dose locality:

~~~text
2.0015
1.0712
1.0091
1.7099
2.2514
1.3188
~~~

Classification:

~~~text
DYNAMIC_NMDA_LOCALITY_ADVANTAGE_PRESENT
~~~

So the quasi-static Gate-15 negative was not the end of the mechanism. When the
authors' released AMPA/NMDA kinetics are restored, compact coactivation on a
single dendritic section gains a robust HUMAN/FROZEN interaction advantage
relative to matched dispersed inputs.

This result does not by itself identify the cause. Gate 17 therefore attacks it
with each exact site's own time-domain superposition null.
