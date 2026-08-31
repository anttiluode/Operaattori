# Gate 6 — out-of-family moving-operator audit

Date: 2026-08-31

Gate 5's positive world was too convenient: a two-mode moving operator was
tested against a hidden two-mode moving teacher. Gate 6 deliberately removes
that exact family match.

The hidden world is now a **six-mode** positive memory kernel. The six time
constants drift around a shared slow trend but also move independently, and
the mixture weights wander among the modes.

The students do not see any teacher time constant or mixture weight.

## Hidden teacher

Base timescales:

~~~text
2.5, 5, 10, 20, 40, 80
~~~

All six are continuously modulated. The global scale follows:

~~~text
0.70 -> 0.95 -> 1.45 -> 1.75 -> 1.05 -> 0.78 -> 1.30 -> 0.72
~~~

and each mode has additional phase-shifted drift.

The positive mixture weights also drift continuously. The target is therefore a
time-varying six-exponential temporal kernel rather than a member of the
two-mode student's exact model family.

## Students and rulers

Moving and frozen diagonal operators are tested at 2, 4 and 8 modes.

A moving N-mode model stores:

~~~text
N hidden states
N eligibility traces
N fast operator states
N slow operator states
N readout weights
----------------------
5N online scalars
~~~

Explicit-context LMS rulers store the recent samples and one adaptive
coefficient per lag.

A deliberately expensive **64-lag recursive least squares** attacker was added
after the first screen so the result could not be dismissed as normalized LMS
adapting too slowly. It keeps a full inverse-covariance matrix: 4,224 online
scalars.

Six seeds, 9,000 steps, first 900 discarded.

## Receipt

| model | online scalars | MSE | operator-horizon correlation |
|---|---:|---:|---:|
| moving-2 | 10 | 3.218e-4 ± 1.5e-5 | **0.889** |
| frozen-2 | 10 | 1.876e-3 ± 8.8e-5 | — |
| **moving-4** | **20** | **3.503e-5 ± 1.5e-6** | **0.860** |
| frozen-4 | 20 | 5.189e-5 ± 1.4e-6 | — |
| moving-8 | 40 | 3.879e-5 ± 1.5e-6 | 0.861 |
| frozen-8 | 40 | 5.354e-5 ± 2.4e-6 | — |
| context-4 | 8 | 2.563e-2 ± 1.8e-3 | — |
| context-8 | 16 | 1.075e-2 ± 5.8e-4 | — |
| context-16 | 32 | 4.308e-3 ± 1.9e-4 | — |
| context-32 | 64 | 1.510e-3 ± 8.5e-5 | — |
| context-64 | 128 | 4.009e-4 ± 3.7e-5 | — |
| context-128 | 256 | 1.942e-4 ± 7.7e-6 | — |
| RLS-64 | 4,224 | 3.959e-4 ± 2.9e-5 | — |

The budget/error Pareto frontier is:

~~~text
context-4   8 scalars    2.563e-2
moving-2   10 scalars    3.218e-4
moving-4   20 scalars    3.503e-5
~~~

Nothing above 20 online scalars improves the mean MSE in this screen.

## What survives

The result from Gate 5 was not merely "a two-mode student rediscovers a
two-mode teacher."

A four-mode moving operator approximates a richer six-mode drifting temporal
law substantially better than its frozen four-mode counterpart and better than
all explicit-context rulers tested here, while keeping only 20 online scalars.

The two-mode moving operator also improves dramatically over its frozen
counterpart.

So the useful statement strengthens to:

> **A low-dimensional operator whose recurrence itself adapts online can act as
> a compact moving approximation to a richer smoothly drifting temporal
> process.**

## Important negative / boundary results

### More moving modes did not help

Eight modes are slightly worse than four:

~~~text
moving-4    3.503e-5
moving-8    3.879e-5
~~~

The update law has a capacity/adaptation tradeoff. "Make the matrix bigger" is
not automatically the answer.

### RLS does not rescue explicit context here

The full-covariance 64-lag RLS attacker is much more expensive and still
performs worse than the 20-scalar moving-4 operator.

That makes the compactness result harder to dismiss as a weak-LMS artifact.

### But this is still a smooth-kernel family

The hidden teacher is richer, but it is still built from sums of exponential
memory modes. Gate 5 already showed an exact-delay world where explicit
attention is decisively superior.

Gate 6 therefore does **not** establish a universal advantage over
Transformer-style context.

The next mathematical weakness is now obvious: all Operaattori moving matrices
so far are diagonal in a fixed basis. Only the eigenvalues move.

A genuinely moving matrix should also be able to move its **basis**.

## Run

~~~bash
python experiments/gate6_rich_operator.py
python -m unittest tests.test_rich -v
~~~

## Stopping line

> **Gate 6 survives model-family mismatch and a strong adaptive-filter
> attacker, but only for smoothly compressible temporal structure. The next
> gate must move eigenvectors, not merely eigenvalues.**
