# Gate 18 — does the effect really follow the human NMDA gamma?

Gate 17 survives a same-site independent-superposition attack, but its
FROZEN comparison is intentionally extreme. Gate 18 replaces that attacker
with the synaptic control already supplied by Aizenbud et al. and their
released FCI code.

## Paper Hybrid B

The released parameter table defines Hybrid B (`human_rat_gamma`) as:

    human AMPA kinetics
    human NMDA kinetics
    human AMPA conductance
    human NMDA conductance
    gamma = 0.062 /mV

while HUMAN differs only by:

    gamma = 0.078 /mV

Thus the first comparison is not a homemade synapse model.

## Harder rest-matched gamma control

Changing gamma also changes the amount of effective NMDA conductance present
at resting voltage. To separate that from the shape of voltage dependence,
Gate 18 adds a harder control:

    gamma = 0.062

but rescales each site's raw NMDA ratio so that at that site's actual settled
pre-event voltage its effective NMDA conductance is exactly equal to HUMAN.

After the event, the two conditions differ because their voltage-dependence
curves have different slopes/shapes, not because one started with more
effective NMDA conductance.

## Ruler

Gate 18 reuses Gate 17's same-site superposition attack at the locked
48-synapse dose.

For each exact three-site set:

    I = AUC(simultaneous response) / AUC(sum of its three single-site traces)

Then compare HUMAN with Hybrid B or the rest-matched gamma=0.062 condition
inside clustered and dispersed arrangements, and finally divide those gains.

The >1.05 threshold is inherited unchanged from Gate 17.

## Outcomes

- HUMAN_GAMMA_LOCALITY_SURVIVES_RESTMATCHED_HYBRID_B
- PAPER_HYBRID_B_DIFFERS_BUT_RESTMATCHED_SLOPE_ATTACK_FAILS
- NO_ROBUST_HUMAN_GAMMA_LOCALITY_ADVANTAGE
- HYBRID_B_ATTACK_MATCH_INADEQUATE

## Stopping line

Only survival of the rest-matched gamma control earns a timing/order
perturbation. It still does not earn growth.


## Receipt — the gamma-slope effect survives, narrowly

The six Gate-17 branches were rerun with two gamma=0.062 attackers.

### Paper Hybrid B

Hybrid B keeps the human synaptic kinetics and raw conductances but substitutes
the smaller rat NMDA gamma.

~~~text
median compact HUMAN / Hybrid-B
  interaction gain                     2.0426x

median dispersed HUMAN / Hybrid-B
  interaction gain                     1.2296x

median locality                        1.6439
fraction locality >1.05                0.833  (5 / 6)
~~~

Per-branch locality:

~~~text
1.8234
1.2610
0.9968
1.7968
2.0757
1.4910
~~~

### Rest-matched gamma=0.062

The harder control rescales the gamma=0.062 raw NMDA ratio separately at every
selected site so its effective NMDA conductance at that site's actual settled
voltage exactly matches HUMAN before the event.

The required raw NMDA ratio is only about one third of the human raw ratio:

~~~text
median ratio scale vs HUMAN            0.3337
~~~

Despite that rest match:

~~~text
median compact HUMAN / restmatched
  interaction gain                     1.5209x

median dispersed HUMAN / restmatched
  interaction gain                     1.2750x

median locality                        1.0592
fraction locality >1.05                0.500  (3 / 6)
~~~

Per-branch rest-matched locality:

~~~text
1.3250
0.8247
0.9605
1.0073
1.2557
1.1110
~~~

Classification:

~~~text
HUMAN_GAMMA_LOCALITY_SURVIVES_RESTMATCHED_HYBRID_B
~~~

This is deliberately recorded as a **narrow pass**, not a universal property.
The paper's Hybrid B separation is large. Once resting effective NMDA strength
is matched away, the median advantage shrinks to about six percent and is
present by the preregistered >1.05 ruler in only three of six branches.

What survives is therefore a specific claim: the shape/steepness of the human
NMDA voltage dependence contributes to compact-branch interaction beyond its
resting effective conductance. The effect is heterogeneous.

## Revised stopping line

This earns one timing/order experiment in the same pinned released model.

Do not introduce growth yet. The next question is whether the spatial order of
the same compact sites and the temporal order of their arrivals interact in a
way that survives each order's own independent-superposition prediction.
