# Operaattori

**Operaattori** is Finnish for **operator**.

This repo asks one deliberately narrow question:

> **Can identical initially matched substrates convert the temporal order of a
> signal into persistent grown anatomy, after all labile state has washed out,
> and does that anatomy causally become a different future operator?**

The headline experiment is the one a wet neural culture cannot easily give us:
**identical-clone order tests with complete internal observability**.

This is not an organoid simulator, a biological-neuron model, or an
intelligence claim. The first object is much smaller:

~~~text
signal history
     ↓
fast excitable state
     ↓
local transport / eligibility
     ↓
slow hysteretic material
     ↓
material changes later transport
     ↓
grown operator
~~~

The matrix is present from commit one. It is the ruler, not the enemy.

## The protocol is designed against the obvious cheat

The primary histories are:

~~~text
H0 = A A A B B B
H1 = A B B A A B
~~~

They have:

- exactly the same number of A and B blocks;
- the same **first** block;
- the same **last** block;
- the same total drive;
- the exact same long suffix after the histories differ;
- a silence interval longer than seven eligibility time constants.

So a decoder cannot win merely by reading "what happened last."

The canonical decoder is deliberately stupid: leave-one-clone-pair-out
**nearest centroid** on the final morphology. A second readout normalizes every
anatomy to unit total mass first, so total material amount cannot be the only
answer.

## The multistability null is load-bearing

Hysteresis is necessary to escape ordinary contraction, but it can also latch
its own noise. Therefore every order result is compared with repeated runs of
the **same history** under independently resampled microscopic noise.

The order effect must clear that within-history multistability floor.

## The rulers

Two same-capacity abstract controls are always reported:

1. **ContractiveMatrixRuler** — 32 stable linear state variables with the same
   two input channels and the same washout. It is allowed to remember only
   through ordinary decaying state.
2. **BistableStateRuler** — 32 abstract hysteretic state variables. This is the
   dangerous null: if hysteresis alone stores the sequence, it should pass.

If the bistable abstract state passes, that is not a failure of the assay. It
prevents the repo from claiming that morphological memory is a computational
function unavailable to ordinary abstract state.

## Gate ladder

### Gate 0 — protocol + ruler

Audit the order protocol before trusting any anatomy. Run both matrix/state
rulers.

### Gate 1 — signal → morphology

After identical suffix + long silence:

- fast state must be gone;
- labile eligibility must be gone;
- raw anatomy must decode history;
- **mass-normalized anatomy** must still decode history.

This earns morphological memory, not learning or intelligence.

### Gate 2 — multistability noise floor

Repeat each history with independently resampled noise and initial
micro-heterogeneity. The between-history effect must exceed the latch's own
within-history spread. A same-history pseudo-classification must stay near
chance.

### Gate 3 — anatomy → operator

Freeze anatomy, reset every fast variable, and apply the same standardized
probe pulses.

The two histories are pairwise mass-matched before probing. If their responses
still differ, spatial arrangement is causal.

For this first substrate the frozen fast dynamics are deliberately linear, so
the repo also derives the exact matrix

~~~text
u(t+1) = A[m] u(t) + B s(t)
~~~

from the grown anatomy and verifies that matrix replay matches the spatial
simulation numerically.

That is the point of the name:

> **the anatomy is a self-grown operator.**

It does **not** yet beat a matrix. When frozen, it *is exactly representable as
one*.

### Gate 4 — LOCKED until Gates 1–3 pass

Only after the first three gates survive: add soma consequence and ask whether
useful morphology is stabilized by delayed local credit.

### Gate 5 — LOCKED behind Gate 4

Only then run functional attackers: proportional-growth null, matched random
arbor, shuffled anatomy, fixed reservoir, learned matrix, and a fixed-capacity
nonlinear surrogate.

Do not build a bug, animal, population, genome, reproduction loop, or
"intelligence" layer before the single substrate earns itself.

## Why this exists

The direct ancestry is spread across earlier repos:

- [FunctionalArbors](https://github.com/anttiluode/FunctionalArbors) — free
  branching material could grow task-relevant geometric delay; later gates
  exposed credit assignment as the hard part.
- [Sunday](https://github.com/anttiluode/Sunday) — history can persist in
  distributed substrate state under a material budget; reservoir attackers
  demoted the broad computing claim.
- [GeoNeuronX](https://github.com/anttiluode/GeoNeuronX) — dendritic dynamics
  materialize temporal history into simultaneous branch state; local
  nonlinearities act before soma collapse.
- [yrotisopeRweN](https://github.com/anttiluode/yrotisopeRweN) and
  [T-800NNP](https://github.com/anttiluode/T-800NNP) — continuously running
  state, stable local addresses, delayed eligibility and finite structural
  allocation.
- [DifferentMachine](https://github.com/anttiluode/DifferentMachine) —
  inherited developmental rule is distinct from acquired structural
  phenotype.
- [BlackBoxLab](https://github.com/anttiluode/BlackBoxLab) Datarium 5 —
  activity writes slow oriented internal material and that material changes
  later field geometry.

Operaattori is an attempt to put the smallest surviving pieces into one
controlled object rather than adding another life simulation around them.

## Run

~~~bash
python -m pip install -r requirements.txt

python experiments/gate0_protocol_and_rulers.py
python experiments/gate1_order_memory.py
python experiments/gate2_multistability_null.py
python experiments/gate3_grown_operator.py

python -m unittest discover -s tests -v
~~~

Machine-readable receipts are written under `results/`.

## First receipt — Gates 0–3 survive the initial implementation

The first 24-clone / 12-clone battery gives:

| measurement | result |
|---|---:|
| contractive matrix order decoder | **0.500** |
| bistable abstract-state order decoder | **1.000** |
| raw morphology order decoder | **1.000** |
| **unit-total-mass morphology** order decoder | **1.000** |
| fast-state residual after washout | **2.43e-7** |
| eligibility residual after washout | **2.04e-3** |
| order-shape distance / same-history multistability floor | **125.67×** |
| same-history pseudo-class decoder | **0.417** |
| mass-matched future-response order distance | **1.028 relative L2** |
| response order effect / same-history response floor | **23.26×** |
| mass-matched frozen-operator matrix distance | **0.0178 relative L2** |
| exact spatial-vs-matrix replay error | **1.05e-17** |

This is a useful combination of positive and negative results.

**Positive:** order survives as spatial morphology after the labile state is
gone, survives total-mass normalization, clears the hysteresis/noise floor, and
changes later standardized propagation after pairwise total-mass matching.

**Negative / ruler:** a same-capacity abstract bistable state also stores the
order perfectly. And once the morphology is frozen, the spatial probe is
exactly reproducible by the derived linear matrix. Operaattori has therefore
earned **self-grown operator**, not a new function class and not superiority to
matrices.

Detailed receipts: [Gate 0](results/GATE0.md),
[Gate 1](results/GATE1.md), [Gate 2](results/GATE2.md), and
[Gate 3](results/GATE3.md).

## Current stopping line

> **Gates 1–3 now earn morphological memory and a causal grown operator in the
> first controlled strip. Gate 4 may now ask whether delayed consequence can
> stabilize useful structure. Gate 5 remains locked: no functional-complexity
> or intelligence comparison is allowed until consequence-guided growth
> survives.**
