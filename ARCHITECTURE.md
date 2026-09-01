# Operaattori architecture

The current reusable hypothesis is a composition, not a new neuron mechanism.

~~~text
input / event pattern
        |
        v
local nonlinear branch operator N_b
        |
        | site-wise transmembrane synaptic currents
        v
geometry-dependent transport operator T_g
        |
        v
somatic response
~~~

The scientific question is whether these two objects can be changed and reused
independently.

The first architecture audit is preregistered in
[results/OPERATOR_FACTORIZATION.md](results/OPERATOR_FACTORIZATION.md).

Its strongest test is deliberately out-of-distribution with respect to the
local operator: N_b is measured only in the original released cell, while T_g
is remeasured after plus/minus 20 percent intrinsic branch-length perturbations.
No held-out somatic trace is used to fit the prediction.

If it works, geometry becomes a replaceable transport module rather than a
reason to relearn the local nonlinear computation.

If it fails, the oracle decomposition tells us whether the failure belongs to
transport linearization or to portability of the local nonlinear operator.
