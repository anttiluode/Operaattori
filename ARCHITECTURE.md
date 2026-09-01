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


## Cross-input reuse

A second audit reused each held-out geometry's single T_g across three distinct
local input operators: one middle site, an outer-site pair, and all three
sites.

~~~text
36 held-out cases

frozen-soma attacker                  0.0931 median NRMSE
factorized T_g[N_original(pattern)]   0.0282
transport oracle                      0.0039

factorized wins                       28 / 36

pattern medians
  middle single                       0.0125
  outer pair                          0.0373
  triple                              0.0512
~~~

Classification:

~~~text
TRANSPORT_OPERATOR_REUSES_ACROSS_INPUT_PATTERNS
~~~

This strengthens the modular claim: T_g is not tied to the waveform used to
measure it.

The remaining limitation is now clearer. N_b is still supplied as a measured
original-geometry current waveform for each input pattern. The next reduction
should replace that lookup with a geometry-independent synaptic nonlinearity
coupled through a local geometry-dependent Green/impulse matrix.


## Reduced nonlinear circuit

The per-pattern local-current lookup can now be removed entirely.

A three-site reduced circuit was built from:

~~~text
released AMPA/NMDA conductance law
        x
3x3 local Green matrix G_g
        x
site-to-soma transport T_g
~~~

The nonlinear site currents are solved self-consistently:

~~~text
V_local = V_noinput + G_g * J
J       = synapse_law(V_local, input)
V_soma  = T_g * J
~~~

Across 54 full-model comparisons spanning six branches, three geometries
(original, 0.80x length, 1.20x length) and three input patterns:

~~~text
transport oracle soma NRMSE      0.0040
reduced circuit soma NRMSE       0.0043
reduced current NRMSE            0.0038

held-out frozen-current model    0.0282
held-out reduced circuit         0.0043
reduced / frozen-current         0.1527
reduced wins                     32 / 36
~~~

Classification:

~~~text
LOCAL_GREEN_MATRIX_X_SYNAPSE_LAW_REDUCES_RELEASED_NEURON
~~~

This is the current reusable object of Operaattori.

The full released branch response no longer needs to be cached as a nonlinear
waveform. Its behavior is reconstructed from an explicit local synapse law
embedded in measured linear transport.

The next attack is temporal rather than geometric: reuse the same reduced object
for asynchronous event patterns without fitting new conductance shapes or
transport kernels.


## Temporal portability

The same reduced circuit was then tested on asynchronous event programs without
new conductance or transport fitting.

~~~text
72 full-model comparisons

transport oracle soma NRMSE       0.0037
reduced circuit soma NRMSE        0.0050
reduced current NRMSE             0.0026

timing medians
  synchronous                     0.0070
  forward 0,5,10 ms               0.0053
  reverse 10,5,0 ms               0.0054
  spread 0,15,30 ms               0.0041
~~~

Classification:

~~~text
TEMPORAL_GREEN_CIRCUIT_GENERALIZES_WITHOUT_REFIT
~~~

A diagnostic open-loop version that keeps the same conductance templates and
transport but removes local voltage feedback has **0.4141 median soma NRMSE**.
The nonlinear fixed-point circuit beats it in **72 / 72** cases.

So the current factorization is genuinely mixed:

~~~text
linear geometry-dependent transport
             x
essential local nonlinear voltage feedback
~~~

The next missing map is upstream of both:

~~~text
geometry parameters  --->  G_g, T_g
~~~

At present those kernels are still measured separately for each geometry.
