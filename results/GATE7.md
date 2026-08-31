# Gate 7 — move the basis

Date: 2026-08-31

Gates 5-6 moved recurrence eigenvalues while keeping the operator basis fixed.
That is only half a moving matrix.

Gate 7 asks whether the **eigenvectors themselves** can become dynamical state.

The hidden teacher is:

~~~text
A(t) = R(phi(t)) diag(lambda_1(t), lambda_2(t)) R(phi(t))^T
~~~

with two drifting timescales and a continuously rotating basis.

The learner sees a two-dimensional input, predicts the hidden teacher's current
two-dimensional state, and receives that target only after prediction.

## Online sensitivity

For operator coordinate theta_j the learner carries the forward sensitivity

~~~text
P_j(t) = d h(t) / d theta_j
~~~

with recurrence

~~~text
P_j(t+1)
  = A(t) P_j(t)
  + (dA/dtheta_j) (h(t) - x(t))
~~~

There is no replay buffer and no BPTT.

The three operator coordinates are:

~~~text
logit(lambda_1)
logit(lambda_2)
phi
~~~

## Gate 7 v1 failed

The first implementation normalized each coordinate's credit independently.

It improved function:

~~~text
moving-full MSE        0.007958
moving-diagonal MSE    0.012389
frozen MSE             0.024586
~~~

and improved direct operator reconstruction:

~~~text
full operator error    0.0338
diagonal error         0.0626
~~~

but it failed the preregistered basis criterion:

~~~text
basis alignment        0.146
required              >= 0.600
~~~

That failure remains in the kill ledger.

## Why v1 failed mathematically

The eigenvalue and angle sensitivities are coupled. Treating them as three
independent scalar credit channels was not justified.

The second attempt kept the same gate thresholds and changed only the local
coordinate solve.

Instead of diagonal normalization, it uses the tiny coupled sensitivity system:

~~~text
(P^T P + 0.02 I) delta = P^T error
~~~

This is a three-coordinate online Gauss-Newton-like step. It still uses only
the current local sensitivity matrix and current consequence; it does not
differentiate through the historical sequence.

## v2 receipt

Six seeds, 9,000 steps, first 900 discarded:

| model | online scalars | MSE | direct operator error | basis alignment |
|---|---:|---:|---:|---:|
| **moving full basis** | **14** | **0.004382 ± 0.00066** | **0.0239** | **0.671** |
| moving diagonal | 10 | 0.012389 ± 0.00065 | 0.0626 | — |
| frozen full | 14 | 0.024586 ± 0.00058 | — | — |
| context-4 | 24 | 0.037858 ± 0.0010 | — | — |
| context-8 | 48 | 0.016389 ± 0.00040 | — | — |
| context-16 | 96 | 0.006627 ± 0.00030 | — | — |
| context-32 | 192 | 0.002315 ± 0.00017 | — | — |
| RLS-32 | 4,288 | 0.001826 ± 0.00011 | — | — |

The unchanged preregistered criteria now pass:

~~~text
full MSE <= 0.80 * diagonal MSE      PASS
full MSE <= 0.80 * frozen MSE        PASS
full operator error <= 0.80 * diag   PASS
basis alignment >= 0.60              PASS (0.671)
~~~

## Budget/error frontier

~~~text
moving-diagonal     10 scalars    0.012389
moving-full         14 scalars    0.004382
context-32         192 scalars    0.002315
RLS-32            4288 scalars    0.001826
~~~

The moving full operator is not the lowest-error model. Explicit context wins
once given substantially more state.

But allowing the **basis itself** to move creates a large gain over the
fixed-basis moving attacker at very small online state.

## What Gate 7 earns

The phrase "moving matrix" is no longer shorthand for moving diagonal decay
constants.

The learner now changes:

- its eigenvalues;
- its eigenvectors;
- therefore which combinations of state are slow/fast directions.

The direct god-mode operator reconstruction improves at the same time as the
functional prediction.

So:

> **The moving operator can alter both temporal scale and computational basis
> online from local forward sensitivity plus delayed consequence.**

## What it does not earn

- a Transformer replacement;
- arbitrary full-matrix plasticity;
- nonlinear computation;
- a learned developmental law;
- superior absolute accuracy to large explicit context.

The update law is still hand-designed.

That is now the obvious next question.

## Run

~~~bash
python experiments/gate7_rotating_basis.py
python -m unittest tests.test_rotating -v
~~~

## Stopping line

> **Gate 7 earns a genuinely moving operator basis. The next rung is no longer
> "make A more dynamic." It is: can one fixed learned developmental law G_theta
> generate useful A(t) trajectories across worlds it was not tuned on?**
