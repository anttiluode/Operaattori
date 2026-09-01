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


## Earned result

The first held-out composition audit passed.

~~~text
original transport reconstruction         0.0046 median NRMSE

held-out +/-20% branch-length perturbations
  frozen original soma attacker           0.0982
  T_g [ N_original ]                      0.0512
  T_g [ N_holdout ] transport oracle      0.0046

factorized / frozen median error           0.5213
factorized beats frozen                    9 / 12
median local-current waveform drift        0.0560
~~~

Classification:

~~~text
TRANSPORT_X_LOCAL_NONLINEAR_OPERATOR_FACTORIZATION
~~~

So the useful abstraction is no longer just a diagram. On these six real
cell-1125 compact branches, a nonlinear synaptic-current operator measured only
in the original geometry can be reused after held-out intrinsic metric changes
by replacing the site-to-soma transport operator.

The transport oracle is especially strong: once given the correct held-out
currents, the linear transport module reconstructs soma traces at about 0.46%
median NRMSE. The dominant residual error in the reusable factorization is
therefore drift of the local nonlinear current operator, not failure of the
transport decomposition.

This does not mean N_b is universal. The next audit asks whether the same T_g
works across several distinct local input patterns rather than only the
three-site Gate-20 pattern.


## Literature fence

The broad architecture is not a novelty claim.

Poirazi, Brannon & Mel (2003), *Pyramidal neuron as two-layer neural network*,
showed that a detailed CA1 pyramidal model could be abstracted as nonlinear
dendritic subunits whose outputs are pooled before somatic thresholding:

https://doi.org/10.1016/S0896-6273(03)00149-1

Subsequent experimental and modeling work developed and attacked this
two-layer / nonlinear-subunit picture, including branch independence and
location-dependent subunit behavior.

Likewise, linear cable transfer and Green's-function descriptions of current
propagation from dendritic sites to the soma are classical cable theory.

Operaattori therefore does **not** claim to invent nonlinear dendritic subunits,
linear transport, or a two-layer neuron.

The narrower earned object here is a causal modularity assay on one released
human L2/3 reconstruction:

~~~text
measure N only in original geometry
change intrinsic branch metric
remeasure only T
predict unseen held-out soma trace with T_new[N_original]
compare against both a frozen-output attacker and T_new[N_actual] oracle
~~~

That intervention-based portability test is what the current receipts support.
Priority beyond that would require a dedicated literature review.
