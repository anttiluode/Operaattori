# Gate 17 — attack Gate 16 with each set's own independent sum

Gate 16 is positive, but passive DC matching cannot guarantee that the
clustered and dispersed sites have identical passive **temporal** filtering.

Gate 17 therefore removes that explanation almost entirely.

For one exact three-site set and one condition:

1. activate all three sites simultaneously;
2. activate site A alone and record voltage at all three sites;
3. activate site B alone and record voltage at all three sites;
4. activate site C alone and record voltage at all three sites;
5. sum the three single-site depolarization traces.

The summed trace is the set's own independent-superposition prediction. It
already contains that set's individual input resistance, cable filtering,
soma transfer, and passive voltage spread.

Define:

    I = AUC(simultaneous local response) / AUC(sum of single-site traces)

and then the NMDA interaction gain:

    G = I_HUMAN / I_FROZEN

Finally compare compact and dispersed sets:

    L_interaction = G_clustered / G_dispersed

The frozen control is tightened relative to Gate 16: every selected synapse is
frozen at its own actual settled pre-event voltage, not a universal -70 mV.

Gate 17 runs only the high 48-synapse dose. It is an attacker, not a new search
for a threshold.

## Stopping line

If the Gate-16 locality effect disappears here, individual temporal filtering
explained it and the claim closes.

If it survives, the candidate mechanism is genuinely interaction-shaped: human
NMDA changes the departure from each exact set's independent superposition more
strongly for compact same-branch inputs. That would earn a timing-structure
experiment, not developmental growth.


## Receipt — survives its own single-site time traces

Gate 17 was expanded from four to all six Gate-16 branches without changing
the metric or threshold.

~~~text
branches                               6
median passive Rinput match            1.037x
median passive soma-transfer match     1.126x
settled selected-site voltage          -70.485 .. -70.400 mV

median compact HUMAN/FROZEN
  interaction gain                     2.6265x

median dispersed HUMAN/FROZEN
  interaction gain                     1.3818x

median interaction locality            1.6759
fraction locality >1.05                0.667  (4 / 6)
median peak-interaction locality       1.0313
~~~

Per-branch locality:

~~~text
1.8765
0.9850
0.9780
1.6952
2.1013
1.6565
~~~

Classification:

~~~text
NMDA_INTERACTION_LOCALITY_SURVIVES_SUPERPOSITION_ATTACK
~~~

This is a stronger statement than Gate 16 but still a heterogeneous one. Two
branches are essentially null. Four show a substantial locality effect.

Because each simultaneous response is divided by the sum of that exact set's
own three single-site time traces, differences in individual passive temporal
filtering are not sufficient to explain the median effect.

Gate 18 therefore replaces the extreme frozen-NMDA attacker with the paper's
own Hybrid B control and a still harder rest-matched gamma=0.062 control.
