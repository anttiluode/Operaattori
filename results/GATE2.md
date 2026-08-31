# Gate 2 — multistability must beat its own noise

## Question

Is Gate 1 merely a bistable medium remembering its own microscopic coin flips?

## Design

For 24 clone seeds:

1. paired H0 vs H1 runs use the same starting clone and the same microscopic
   noise schedule;
2. the same history is then rerun with independent microscopic noise to
   measure the latch's own path-dependent spread;
3. two independently noisy H0 repeats are falsely labeled as separate classes
   and given to the same nearest-centroid decoder.

All distances use unit-total-mass anatomy.

## Receipt

~~~text
mean paired H0/H1 shape distance      0.048682
mean same-history latch/noise floor   0.000387
order / multistability floor            125.67 x

same-history pseudo-class accuracy       0.417
~~~

## Result

The signal-order effect is far larger than the morphology variation produced
by rerunning the same history through independently noisy copies. The stupid
decoder cannot reliably classify which noise stream produced an otherwise
identical history.

## Stop line

> **The bistable medium is path dependent, but its self-generated path
> dependence is not large enough to explain the order separation.**
