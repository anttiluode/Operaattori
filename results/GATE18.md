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
