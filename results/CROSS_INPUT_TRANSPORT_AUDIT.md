# Cross-input transport reuse audit

The first operator-factorization audit passed on the Gate-20 three-site input.
This attacker asks whether the measured transport operator is genuinely
input-independent.

## Locked object

Reuse:

- the same six Gate-20 branches;
- the same HUMAN AMPA/NMDA mechanisms;
- multiplicity 8 per active site;
- held-out branch-length scales 0.80 and 1.20;
- the same matched-control 0.001 nA one-dt transport kernels.

For each branch define three input patterns by the ordered compact sites:

~~~text
middle_single   = [middle]
outer_pair      = [proximal, distal]
triple          = [proximal, middle, distal]
~~~

No pattern is selected after seeing performance.

## Factorization

For every pattern p, measure its local synaptic-current output only in the
original geometry:

~~~text
N_b(p) -> J_original,p(t)
~~~

For each held-out geometry g, measure the three site-to-soma transport kernels
once, forming one T_g.

Then use that same T_g for every pattern:

~~~text
Vhat_g,p = T_g [ J_original,p ].
~~~

The held-out nonlinear current J_g,p is recorded only for the diagnostic
transport oracle.

## Cases

~~~text
6 branches
x 2 held-out metric scales
x 3 input patterns
= 36 held-out predictions
~~~

## Locked rulers

The audit earns:

~~~text
TRANSPORT_OPERATOR_REUSES_ACROSS_INPUT_PATTERNS
~~~

only if:

1. median original-geometry reconstruction NRMSE <= 0.10;
2. median held-out transport-oracle NRMSE <= 0.10;
3. median factorized held-out NRMSE <= 0.15;
4. factorized median NRMSE <= 0.80 x frozen-soma attacker median NRMSE;
5. factorized beats frozen in at least 24 / 36 holdouts;
6. each of the three pattern families has median factorized NRMSE <= 0.15;
7. no actual held-out run reaches the -20 mV soma spike guard.

No trace alignment, per-pattern kernel, fitted gain, branch selection, or input
dose change is allowed.

Passing means T_g is a reusable geometry module across distinct local operator
outputs. It does not yet mean N_b can predict an input pattern it has never
seen.

This is an architecture audit, not Gate 25.
