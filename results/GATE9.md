# Gate 9 — history as a path through operator space

The trigger was the matrix-exponential observation:

```text
dx/dt = A x
x(t)  = exp(A t) x(0)
```

For a moving generator, the correct object is a time-ordered product. If two
successive generators do not commute, order can matter even when their total
first-order exposure is identical.

Gate 9 audits that mechanism directly before attaching another organism,
tissue, or learning story to it.

## Primitive generators

Only two 2×2 nilpotent shears are supplied:

```text
H = [[0, 1],      V = [[0, 0],
     [0, 0]]           [1, 0]]
```

Their commutator is

```text
[H,V] = H V - V H
      = [[ 1, 0],
         [ 0,-1]]
```

—a saddle direction that was **not** a primitive.

## Closed-loop assay

Run the operator path

```text
H -> V -> -H -> -V
```

and compare it with the reversed loop. The integrated generator is exactly
zero, so an order-blind approximation based only on

```text
exp(integral A dt)
```

predicts identity.

Measured at epsilon = 0.08:

| measurement | result |
|---|---:|
| net integrated-generator norm | **0** |
| noncommuting closed-loop residue | **0.0091088** |
| commuting-control residue | **3.14e-16** |
| reverse-loop difference alignment with `[H,V]` | **0.996835** |
| residue scaling exponent vs epsilon | **2.0029** |

The epsilon-squared scaling is the expected Lie-bracket / BCH signature. The
path leaves a second-order residue even though the first-order integral is
zero.

## Held-out path-length assay

A sequence is built from clockwise and counter-clockwise closed loops. Every
individual loop uses the same primitive operators and has zero net generator
exposure; only its orientation differs.

A two-feature linear observer sees only the final 2D state.

Training loop counts:

```text
8, 12, 16
```

Held-out loop counts:

```text
10, 14, 18
```

Across six seeds:

| measurement | result |
|---|---:|
| final-state correlation with signed path area | **0.9717 ± 0.0007** |
| held-out R² | **0.9130 ± 0.0071** |
| commuting-control R² | **-0.0006 ± 0.0004** |

So the present state can carry a compact trace of the operator path even after
the instantaneous generator has returned to neutral.

The killer stays beside it: a tiny digital counter computes this toy signed-loop
target exactly. Gate 9 is therefore a **mechanism result**, not a state-efficiency
or superiority claim.

## Lie-closure audit: does a layer above actually appear?

The 2D answer is both positive and limiting.

```text
primitive span:  H, V                 dimension 2
Lie closure:     H, V, [H,V]          dimension 3
```

Then it saturates. Fixed 2×2 linear dynamics cannot keep birthing new algebraic
directions forever.

The more interesting audit uses a three-state chain with only four
nearest-neighbor bidirectional shear primitives:

```text
1 <-> 2 <-> 3
```

Their Lie closure is:

```text
primitive dimension: 4
closure dimension:   8
full sl(3):           8
```

The non-neighbor direction `E13` has normalized distance **1.0** from the
primitive span but **4.44e-16** from the generated closure.

That is the clean version of the "layer above" intuition:

> ordered compositions of local operators can generate effective directions
> that no local primitive directly supplied.

It is still not a new function class relative to an unconstrained full 3×3
matrix; the full matrix can represent those directions directly. The possible
advantage is factorization, locality, development, and online reconfiguration.

## Stopping line

Gate 9 earns:

> **history can live in the geometry of a time-ordered operator trajectory, not
> only in the current operator or in an explicit context buffer.**

It also supplies the next constraint. A single tiny matrix has a finite Lie
closure. If this is going to become the spatial / dendritic object we have been
circling, the next experiment should distribute small local generators over a
graph or field and ask whether **changing route/topology plus local operator
state** creates useful effective operators under a strict local budget.
