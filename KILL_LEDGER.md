# Operaattori kill ledger

This file exists so attractive old interpretations do not quietly return.

## Already killed or demoted in the ancestry

### "More fractal = more intelligent"

No. Box-count slopes are not a functional-complexity objective. Scale span,
surface area, nonlinear compartmentalization and actual input/output function
must be measured separately.

### "A pretty growing arbor is already a computer"

No. DendriticAttentionSystem grew signal-correlated structure, but the geometry
was not sufficiently load-bearing in later computation.

### "Geometry beats digital temporal features"

No. GeoNeuronX repeatedly found explicit delays / FIR / digital RC banks that
matched or beat the morphology on the toy temporal tasks. Geometry may
materialize useful history into space; that is a bias, not automatic
superiority.

### "Persistent medium = unique memory primitive"

No. Sunday's fixed-reservoir/readout audit reproduced much of the controlled
relation-computation story. The surviving distinction was *where* the learned
history lived: in the medium rather than only in a trained observer.

### "Ephaptic phase is the neuronal growth teacher"

Not earned. FunctionalArbors found field geometry could help exploration, but
soma-locked absolute phase did not survive the strong controls.

### "Returning reward solves structural credit"

No. FunctionalArbors v0.8/v0.9 transported consequence back successfully, but
free causal assignment remained the bottleneck.

## Operaattori-specific traps

### Recency

No order result counts unless histories have matched marginals, the same first
and last symbol, an exact common suffix, and a washout several labile-state
time constants long.

### Hysteresis self-memory

No morphology result counts unless between-history separation exceeds repeated
same-history multistability/noise spread.

### Total-mass cheat

The canonical order decoder also receives unit-total-mass anatomy. Gate 3
pairwise mass-matches anatomies before probing.

### Matrix avoidance

Forbidden. The matrix is a ruler from Gate 0. Gate 3 explicitly derives the
exact frozen operator matrix and requires it to replay the spatial simulation.

### Gate inflation

Gate 4 did not enter CI until Gates 1–3 survived. Gate 5 was then unlocked only
to test the moving-matrix hypothesis against explicit-context rulers. Gate 6
does not exist yet. No bug, animal, tissue, reproduction, population, language,
reward hierarchy, lesion story, or "intelligence" layer is allowed to jump
this queue.

### "Moving matrix replaces attention"

Killed in Gate 5 for the current two-mode formulation.

A moving recurrent operator tracks a smoothly drifting memory kernel compactly,
but on an exact drifting-delay task causal explicit-token attention wins by
more than an order of magnitude in MSE and tracks the useful lag almost
perfectly. Exact addressable history cannot be assumed to survive compression
into operator state.

The next architecture must therefore treat "compress into current operator" and
"keep explicitly addressable" as distinct memory choices rather than declaring
one universally superior.


### Gate 7 v1 — "moving basis" was not yet identified

The first rotating-basis attempt **failed its preregistered basis criterion**.

Screen:

~~~text
moving-full MSE             0.007958
moving-diagonal MSE         0.012389
frozen-full MSE             0.024586

full operator error         0.0338
diagonal operator error     0.0626

basis alignment             0.146
required                    >= 0.600
~~~

Allowing an angle coordinate clearly improved prediction and direct operator
reconstruction, but the learned basis angle did not track the hidden basis
well enough. This is not counted as a Gate-7 pass.

The likely mathematical defect is that the first update treated the three
operator coordinates independently even though their forward sensitivities are
coupled. The next attempt keeps the threshold and replaces the diagonal
per-coordinate credit normalization with a tiny coupled sensitivity solve.
