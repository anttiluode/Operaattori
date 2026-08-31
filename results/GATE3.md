# Gate 3 — anatomy becomes the operator

## Question

Does the remembered spatial anatomy causally change future signal transport
after every fast variable is reset?

## Design

12 matched clone pairs.

For each H0/H1 pair:

1. develop and wash out as in Gates 1–2;
2. scale the two anatomies to **exactly equal total material**;
3. freeze morphology;
4. reset fast state to zero;
5. apply identical standardized pulses from port A and port B;
6. read the same fixed soma location.

A same-history independently-noisy growth pair supplies the functional noise
floor.

## Receipt

~~~text
mass-matched H0/H1 response distance      1.027988 relative L2
same-history response floor               0.044187
order response / noise floor                23.26 x

mass-matched frozen A[m] matrix distance  0.017800 relative L2

erase spatial arrangement by replacing
both with the same uniform equal-mass
material:
order response distance                   0.000000

maximum direct spatial vs exact matrix
replay error                            1.05e-17
~~~

The response metric is a concatenation of the soma traces produced by the two
standardized probe ports.

## The "organic matrix" result

For frozen material, the fast probe equation is linear:

~~~text
u(t+1) = A[m] u(t) + B s(t)
~~~

`A[m]` is derived directly from the material-dependent local conductances.
The derived matrix reproduces the explicit spatial simulation to numerical
precision.

So the first Operaattori is not a mysterious alternative to a matrix.

The interesting direction is the reverse:

~~~text
signal history
     ↓
local developmental physics
     ↓
morphology m
     ↓
operator matrix A[m]
~~~

The operator was not fitted by backpropagation or supplied as a target matrix.
It was **grown** by the signal/substrate loop.

## Stop line

> **Gate 3 earns a causal self-grown operator. It explicitly does not earn a
> function class that a matrix cannot represent: once frozen, the current
> substrate is exactly linear. Gate 4 must now ask whether consequence-guided
> development can grow useful operators rather than merely different ones.**
