# Gate 4 — delayed consequence shapes the operator

## Question

After Gates 1–3 established persistent morphology and a causal grown operator,
can a delayed scalar soma consequence stabilize a more useful operator through
a lingering local eligibility field?

## Task

The same continuous strip receives a balanced randomized sequence of A and B
trials.

The intentionally simple target is:

~~~text
A input -> larger soma response
B input -> smaller soma response
~~~

The target error is measured at the soma, then withheld for a silent delay.
When it arrives, it can modify material only through the **currently local
eligibility field**. The fast field and eligibility are never hard-reset
between trials.

This task is trivial for an explicit two-weight digital model. The digital
attacker has zero target error by construction. Gate 4 is a credit/growth
test, not a hard-computation benchmark.

## Attackers

~~~text
CAUSAL
correct delayed scalar × correct local eligibility

CREDIT SHUFFLE
same consequence magnitude, randomized sign × correct eligibility

ELIGIBILITY SHUFFLE
correct delayed scalar × spatially permuted eligibility

NO CREDIT
same traffic and material decay, no structural consequence
~~~

12-seed screen, 120 balanced trials:

~~~text
mean selectivity

causal                  +0.481 ± 0.115
credit shuffle          -0.226 ± 0.397
eligibility shuffle     -0.210 ± 0.266
no credit               -0.390 ± 0.055

paired arm-mean differences

causal - credit shuffle        +0.707
causal - eligibility shuffle   +0.692
causal - no credit             +0.871
~~~

## What this earns

The delayed scalar is not sufficient by itself: randomizing its sign destroys
the effect.

The eligibility field is not decorative: keeping the scalar correct while
scrambling the spatial address also destroys the effect.

So in this toy strip:

> **recent local transport leaves an addressable eligibility trace, and a
> later scalar soma consequence can use that trace to change persistent
> material so the final operator becomes more selective.**

## What this does not earn

- hard or general credit assignment;
- biological retrograde signalling;
- superiority to ordinary learning;
- a task that requires morphology;
- a nonlinear neuron;
- intelligence.

The explicit digital solution is two weights. It wins on engineering
simplicity.

The next gate must therefore attack whether the *developmental substrate*
buys something ordinary fixed computation does not: continual regrowth,
damage recovery, partition inheritance, or a function produced without
backpropagating through the machine.

## Stop line

> **Gate 4 earns a first delayed local structural-credit effect. Gate 5 remains
> locked until this survives a broader confirmation and stronger ordinary
> attacker battery.**
