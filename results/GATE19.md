# Gate 19 — does the dendritic scaffold care which way time runs across it?

Gate 9 showed mathematically that order through noncommuting operators can
matter. Gates 16–18 now provide a real reconstructed-neuron substrate with a
measurable compact nonlinear interaction. Gate 19 asks for the literal bridge:

> If the same three compact dendritic sites receive the same events with the
> same amplitudes and delays, does reversing proximal-to-distal versus
> distal-to-proximal arrival order change the nonlinear interaction?

## Locked timing

Three compact sites are sorted by NEURON section x, from proximal to distal.

Primary event times are fixed before the run:

    proximal -> distal: 60, 64, 68 ms
    distal -> proximal: 68, 64, 60 ms

The 4-ms inter-site lag is locked. It is close to the released human NMDA
5-ms rise time but is not searched or optimized.

Each site receives the locked Gate-17/18 high dose of 16 virtual synapses.

## The linear timing attacker

A reversed temporal sequence trivially changes a passive waveform if different
sites have different filters. Gate 19 removes that explanation separately for
every order.

For each exact site set, condition and order:

1. run all three timed events together;
2. run each site alone at its exact event time;
3. sum those three single-site traces;
4. divide simultaneous local-voltage AUC by the timed independent sum.

Thus the forward and reverse interaction ratios are each compared against their
own time-shifted linear prediction.

Order magnitude is:

    | log(I_proximal-first / I_distal-first) |

so a branch is allowed to prefer either direction.

## Controls

Primary condition:

    HUMAN gamma = 0.078

Attacker:

    gamma = 0.062, rest-matched site by site to HUMAN

The same mapped dispersed controls from Gates 16–18 receive the corresponding
event times. That asks whether any order effect is specifically amplified by
compact same-branch geometry.

## Preregistered ruler

The five-percent ruler is reused:

    order magnitude >= log(1.05)

A compact HUMAN order effect must clear that on at least four of six branches.
To call it compact-specific, the compact-minus-dispersed median order magnitude
must also exceed log(1.05), with compact larger on at least four of six
branches.

A stronger final class requires the HUMAN-vs-restmatched gamma contribution to
clear the same aggregate ruler.

## Stopping line

A positive compact-specific result earns a multi-pulse sequence-discrimination
test. It still does not earn developmental growth.
