# Gate 5 — the matrix itself is state

Date: 2026-08-31

Gate 5 is the first explicit test of the "moving matrix" idea.

The hypothesis is not that matrices disappear. It is the opposite:

> **the operator matrix itself becomes part of the dynamical state.**

The minimal model is

~~~text
h(t+1) = A(t) h(t) + (I - A(t)) 1 x(t)
y_hat  = w(t)^T h(t)
~~~

with

~~~text
A(t) = diag(lambda_1(t), lambda_2(t))
logit(lambda) = base + fast(t) + slow(t)
~~~

Prediction error is revealed only after prediction. A local sensitivity trace

~~~text
p_i(t) = d h_i(t) / d logit(lambda_i(t))
~~~

acts as eligibility. The operator update uses that local trace; there is no
BPTT and the moving model keeps no explicit token history.

This gate has **two worlds on purpose**. A positive result counts only beside a
world designed to kill a universal "moving matrices replace attention" claim.

## Online-state budget

The moving model has two modes.

Its online state is:

~~~text
2 hidden states
2 eligibility traces
2 fast operator states
2 slow operator states
2 readout weights
--------------------------
10 mutable/stored scalars
~~~

The fixed base timescales are not counted as online state.

The explicit-context rulers keep both recent tokens and adaptive position
weights:

~~~text
64-token smooth-kernel ruler: 128 online scalars
32-token delay rulers:          64 online scalars
~~~

This is not a parameter-matched benchmark. It is a **history-storage audit**:
can a tiny evolving operator compress a drifting temporal law into its present
state, and where does explicit context remain decisively better?

## World A — smoothly drifting memory kernel

The teacher itself is a hidden two-mode moving operator.

Its two time constants drift continuously and reverse:

~~~text
fast:  2.5 -> 4 -> 7 -> 12 -> 7 -> 4 -> 2.5
slow:    9 ->14 ->24 -> 40 ->24 ->14 -> 9
~~~

The target is

~~~text
y(t) = 0.65 z_fast(t) + 0.35 z_slow(t)
~~~

where each z is an exponential memory mode. The teacher time constants are
never shown to the learners.

All models receive x(t), make a prediction, and only then receive y(t).

### Rulers

1. **MOVING OPERATOR** — two recurrent eigenvalues change online through
   fast/slow operator state.
2. **FROZEN OPERATOR** — identical two-mode state and identical online readout
   update, but the recurrent eigenvalues cannot move.
3. **EXPLICIT CONTEXT LMS** — stores 64 recent scalar tokens and directly
   adapts 64 context weights. On this scalar task this is the linear-attention /
   FIR limit: keep the past, change how it is weighted.

Eight seeds, 9,000 steps, first 800 discarded:

| model | MSE |
|---|---:|
| **moving two-mode operator** | **5.46e-05 ± 3.27e-06** |
| frozen two-mode operator | 1.70e-04 ± 6.76e-06 |
| 64-token explicit context LMS | 6.17e-05 ± 3.25e-06 |

The moving operator uses about **32%** of the frozen model's MSE and about
**88.5%** of the explicit-context model's MSE.

More importantly, its observer-measured effective time constant follows the
hidden teacher's moving temporal scale:

~~~text
corr(operator effective tau, teacher effective tau)
= 0.9921 ± 0.0012
~~~

So in this matched smooth-kernel family, the sentence

> **the matrix moved with the world**

is literal and measurable.

The result is also compact: ten online scalars versus 128 stored/mutable
context scalars for the 64-token ruler.

## World B — exact drifting-delay attack

Now change only the temporal law.

The target is one exact old sample:

~~~text
y(t) = x(t - L(t))
~~~

and the useful lag drifts and reverses:

~~~text
4 -> 6 -> 9 -> 14 -> 20 -> 14 -> 9 -> 6 -> 4
~~~

This is where explicit access to individual old tokens should win.

### Additional transformer-like ruler

**CAUSAL LAG ATTENTION** stores the recent tokens verbatim. After each target is
revealed, every lag position updates an exponentially weighted correlation
score. Prediction softmax-attends over lag positions.

It is deliberately simple, but structurally transformer-like in the relevant
sense: **history remains explicit and the current computation selects from it**.

Eight seeds:

| model | MSE |
|---|---:|
| moving two-mode operator | 1.087 ± 0.015 |
| frozen two-mode operator | 0.997 ± 0.010 |
| 32-token explicit context LMS | 0.253 ± 0.0065 |
| **32-token causal lag attention** | **0.0957 ± 0.0035** |

The lag-tracking audit is even cleaner:

~~~text
corr(moving operator effective tau, true lag)
= -0.062 ± 0.131

corr(attention estimated lag, true lag)
= 0.9953 ± 0.0018
~~~

The moving operator does not merely lose. On this task, moving its two smooth
memory eigenvalues is the **wrong representation**.

The explicit-token system keeps the thing the task actually asks for: a
particular old sample.

## What Gate 5 earns

It earns a narrow mathematical statement:

> **A recurrent operator whose eigenvalues are themselves fast/slow state can
> track a smoothly drifting temporal kernel online and slightly outperform a
> much larger explicit-context linear ruler in that matched regime, while
> storing an order of magnitude fewer online scalars.**

It also earns a much more important boundary:

> **That does not generalize to exact addressable history. When the task asks
> for one specific old token, explicit context/attention wins by more than an
> order of magnitude in MSE and tracks the useful lag almost perfectly.**

## Why this matters for the Transformer question

The useful distinction is no longer "matrix versus no matrix."

It is:

~~~text
explicit-context machine:
past remains available as data

moving-operator machine:
past is compressed into present computation
~~~

Gate 5 shows one regime where the second representation is efficient, and one
regime where that compression destroys exactly the information the task needs.

That is already enough to kill a naive replacement claim.

A serious Transformer alternative would need a **hybrid** or richer moving
operator that can decide what to compress into operator state and what must
remain explicitly addressable.

## Run

~~~bash
python experiments/gate5_moving_operator.py
python -m unittest tests.test_moving -v
~~~

Machine-readable summary and downsampled seed-0 trajectories are written to:

~~~text
results/gate5.json
~~~

## Stopping line

> **Operaattori now has moving-matrix mathematics that genuinely tracks a
> drifting operator family. It is not a Transformer replacement. Exact-delay
> memory is a decisive counterexample for the current two-mode formulation.**
