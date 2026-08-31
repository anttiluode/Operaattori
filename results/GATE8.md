# Gate 8 — one developmental law across unseen worlds

Gate 8 moved the object of selection one level upward.

The outer loop does **not** select a final recurrent matrix. It selects one
small developmental rule

```text
theta = {
  fast operator plasticity rate,
  slow consolidation rate,
  fast-state retention,
  readout adaptation rate
}
```

on eight training worlds. The winning `theta` is then frozen. Every held-out
world starts from the same operator, hidden state, and readout.

The final Gate-8 CI run passed.

| measurement | result |
|---|---:|
| candidate developmental laws | 40 |
| train worlds | 8 |
| held-out worlds | 8 |
| selected candidate | 26 |
| selected train MSE | **1.612e-05** |
| selected held-out MSE | **1.632e-05 ± 8.9e-06** |
| hand-written Gate-6 rule | 3.015e-05 |
| frozen operator | 3.785e-05 |
| median candidate rule | 3.788e-05 |
| cheating per-world oracle | 1.609e-05 |
| selected / oracle | **1.014×** |
| selected / hand rule | **0.541×** |
| selected / median candidate | **0.431×** |
| held-out candidate rank | **1 / 40** |
| held-out horizon correlation | **0.717** |
| std of final operator tau across held-out worlds | **4.920** |

The selected rule was:

```text
fast_rate      = 5.612778932110066
slow_rate      = 0.09648220974900232
fast_retention = 0.9960976679664698
readout_rate   = 0.2182081920788579
```

## What this earns

One fixed developmental law can start from the same `A(0)` and produce useful,
different operator trajectories in previously unseen members of this world
family.

That supports the clean distinction:

```text
genotype-like object   ~= G_theta
developed phenotype    ~= A(t)
```

The final matrix is not what is inherited.

## What it does not earn

The world family is still narrow: smooth six-mode drifting memory kernels.
The outer population is only 40 hand-parameterized candidate rules. This is not
a universal learning law and not a biological genetics claim.

The 64-token LMS ruler is poor on this particular smooth family, which Gate 5
already showed is exactly where moving operators are naturally strong.
Exact addressable-delay tasks remain the explicit-memory boundary.

Gate 8 therefore unlocks the next mathematical question: whether the **path**
through operator space itself carries function that is lost if only the
instantaneous matrix or its time-average is retained.
