# Gate 1 — signal history becomes persistent morphology

## Question

After all fast and eligibility state is washed out, can the temporal order be
read from the grown anatomy with a deliberately stupid decoder?

## Battery

24 clone pairs.

Within each pair:

- identical initial microscopic anatomy;
- identical microscopic noise stream;
- only the history order differs.

The decoder is leave-one-clone-pair-out nearest centroid.

## Receipt

~~~text
raw anatomy decoder accuracy           1.000
unit-total-mass shape decoder           1.000

maximum fast-state residual          2.43e-7
maximum eligibility residual         2.04e-3

mean total material H0               19.521
mean total material H1               18.594

contractive matrix ruler              0.500
bistable abstract-state ruler         1.000
~~~

The total amount of material is not perfectly matched in this first
developmental regime, so raw mass is not a clean order variable by itself.
That is why the canonical shape attacker normalizes each final anatomy to unit
total mass before decoding. It still classifies every held-out pair correctly.

## Result

The order is stored in persistent **spatial arrangement**, not only in
undecayed fast activity or total material amount.

## Stop line

> **Gate 1 earns morphological memory. It does not earn unique computation:
> the bistable abstract ruler remembers just as well.**
