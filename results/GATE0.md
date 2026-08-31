# Gate 0 — protocol and rulers

## Question

Can the order assay itself be made difficult to cheat, and what do matched
abstract controls do under the exact same protocol?

## Histories

~~~text
H0 = A A A B B B
H1 = A B B A A B
~~~

The histories have the same length, the same number of A and B blocks, the
same first symbol, and the same last symbol.

After the histories differ, both arms receive the exact same 120-step
alternating A/B suffix followed by 800 silent steps.

With `dt=0.08` and eligibility decay `0.12`, the silence is **7.68
eligibility time constants**.

## Rulers

24 matched clone pairs, leave-one-clone-pair-out nearest-centroid decoder:

~~~text
contractive linear matrix state   0.500
bistable abstract state           1.000
~~~

The matrix parameters are fixed across clones. Only microscopic initial state
and process noise are resampled.

## Result

The protocol is strong enough to erase ordinary contractive recency, but a
same-capacity abstract hysteretic state stores the order perfectly.

That second result is intentionally dangerous. It means later morphological
memory cannot be sold as a function unavailable to an ordinary abstract
bistable state.

## Stop line

> **The order protocol is not a recency test. Hysteresis itself is already a
> sufficient abstract memory mechanism, so morphology must earn spatial or
> developmental value beyond merely remembering.**
