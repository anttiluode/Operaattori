# Gate 12 — is the biological cable ordering itself unusual?

Gate 11 established a causal fact:

> reordering the same heterogeneous passive cable segments changes transfer.

That is not yet evidence that the biological ordering is special. A dendrite
normally tapers, and moving thick cable from the proximal end to the distal
end can itself produce a large effect.

Gate 12 attacks that interpretation before adding growth, ephaptic coupling,
or learning.

## Question

For the same human L2/3 morphology 1125:

> **Does the biological proximal-to-distal order occupy an unusual part of the
> transfer-function distribution after progressively stronger taper-preserving
> nulls?**

"Unusual" is deliberately not "better." No biological utility objective is
assigned.

## Transfer signature

At 1, 5, 15, 40, 100, and 300 Hz, each path is summarized by:

- sealed distal voltage gain in dB;
- log input-impedance magnitude;
- unwrapped gain phase.

The biological signature is standardized against each permutation ensemble.
Its RMS z-distance from the ensemble mean is compared with the same statistic
for the null members, yielding an empirical upper-tail probability.

## Nulls / attackers

### 1. Full permutation

Shuffle the exact same (length, radius) segment multiset arbitrarily.

### 2. Ideal monotonic taper

Sort the exact same segments thick-to-thin and thin-to-thick. These are strong
deterministic rulers for how much an idealized taper can move the transfer.

### 3. Endpoint-preserving shuffle

Keep the proximal 10% and distal 10% of physical path length fixed. Shuffle
only the middle. If the biological order becomes ordinary here, Gate 11 was
dominated by boundary placement.

### 4. Coarse-profile-preserving shuffles

Shuffle only within fixed proximal-to-distal distance windows of 10, 25, 50,
and 100 um. Segments never cross a window boundary. The coarse radius-vs-
distance profile is retained while fine sequence is randomized.

## Important statistical limitation

Many root-to-tip paths share proximal branches. Path-level tail fractions are
descriptive assays of one reconstructed cell, not independent biological
replicates and not population p-values.

## Outcome classes

The code reports, rather than forces, one of four outcomes:

- REAL_ORDER_NOT_GLOBALLY_UNUSUAL
- ENDPOINT_POSITION_EXPLAINS_MUCH
- COARSE_TAPER_EXPLAINS_MUCH
- FINE_ORDER_REMAINS_UNUSUAL

A scientifically negative outcome still passes CI. CI assertions cover only
the numerical protocol.

## Run

Moderate default:

    python experiments/gate12_biological_order_audit.py

Stronger local run:

    python experiments/gate12_biological_order_audit.py --paths 64 --permutations 1000 --constrained-permutations 300

Outputs:

    results/gate12/gate12.json
    results/gate12/biological_order_nulls.png

## Stopping line

Gate 12 distinguishes **order matters** from **the real order is unusual**.

Even FINE_ORDER_REMAINS_UNUSUAL would not prove optimality, adaptation, or
intelligence. It would only justify asking whether geometry-changing local
rules can move the operator scaffold toward reproducibly different functional
regimes.
