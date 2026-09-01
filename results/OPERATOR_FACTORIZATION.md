# Operator factorization audit

Operaattori now has three empirical separations:

1. rigid re-embedding is an intrinsic cable symmetry;
2. intrinsic metric changes passive transport strongly;
3. the normalized compact-branch nonlinear interaction law is comparatively
   stable under the same 20 percent metric intervention.

This audit asks whether those observations compose into an actually useful
architecture rather than remaining three separate facts.

## Hypothesis

For a compact branch b in geometry g:

~~~text
input pattern
    |
    v
N_b : local nonlinear synaptic operator
    |
    |  site-wise AMPA + NMDA current waveforms
    v
T_g : geometry-dependent linear transport operator
    |
    v
somatic voltage
~~~

The strong form tested here is:

~~~text
V_soma(g, input) ~= T_g [ N_b(g = original, input) ].
~~~

That is: measure the nonlinear local operator once in the original geometry,
then reuse it after an unseen metric perturbation while updating only the
transport operator.

## Why synaptic current is the interface

Local voltage is not a conserved source variable. Synaptic transmembrane
current is.

The released cell-1125 dendrites are passive apart from their point synapses.
A site-specific small-current impulse therefore gives a direct time-domain
transport kernel from that material site to the soma. Summing the convolution
of each synaptic current with its corresponding kernel is the natural
composition.

No regression coefficient is fitted to held-out soma traces.

## Locked branches and input

Reuse the exact six Gate-20 compact branches.

At every branch:

- 3 compact synaptic sites;
- HUMAN AMPA/NMDA kinetics;
- multiplicity 8 per site;
- all three sites activated simultaneously;
- no field;
- no XYZ re-embedding.

The local operator output is the three site-wise waveforms:

~~~text
J_b(t) = -(i_AMPA(t) + i_NMDA(t))
~~~

measured in the original geometry only.

The minus sign is the NEURON current convention: inward excitatory point-process
current is negative, whereas positive IClamp current is depolarizing.

## Held-out geometries

The nonlinear operator is frozen after the original run.

For each branch independently build fresh released models with only that
branch's intrinsic section length changed to:

~~~text
0.80 x original
1.20 x original
~~~

Diameter, topology, normalized synapse x locations, mechanisms and dose remain
unchanged.

The plus-20-percent arm repeats the already-audited scale. The minus-20-percent
arm is a new symmetric holdout. Neither arm is used to fit the nonlinear
operator.

## Transport operator

For each geometry and each of the three sites, inject a one-time-step
subthreshold IClamp pulse:

~~~text
amplitude = +0.001 nA
duration  = dt = 0.025 ms
~~~

Record the somatic depolarization and divide by 0.001. The resulting discrete
kernel h_g,s[k] is the response to a +1 nA one-sample current input at site s.

The predicted soma trace is:

~~~text
Vhat_g[k] = sum_s conv( J_original,s[k], h_g,s[k] ).
~~~

This uses no held-out nonlinear trace.

## Three predictions / attackers

For each of 12 held-out branch x geometry cases report:

A. Frozen-soma attacker

~~~text
Vhat_frozen = actual original-geometry soma trace
~~~

B. Factorized prediction

~~~text
Vhat_factor = T_g [ J_original ]
~~~

C. Transport oracle

Measure the actual held-out synaptic-current waveforms J_g only for diagnosis:

~~~text
Vhat_oracle = T_g [ J_g ].
~~~

The oracle is not a prediction system. It says how much error belongs to the
linear transport approximation versus drift of the local nonlinear current
operator.

## Locked metrics

Trace error:

~~~text
NRMSE = rms(pred - actual) / rms(actual)
~~~

Also report peak and positive-AUC relative error.

Before interpreting held-out results, the same transport composition must
reconstruct the original geometry:

~~~text
T_original[J_original]
~~~

## Locked rulers

The architecture earns:

~~~text
TRANSPORT_X_LOCAL_NONLINEAR_OPERATOR_FACTORIZATION
~~~

only if all of the following hold:

1. median original-geometry reconstruction NRMSE <= 0.10;
2. median held-out transport-oracle NRMSE <= 0.10;
3. median held-out factorized NRMSE <= 0.15;
4. factorized median NRMSE <= 0.80 x frozen-soma attacker median NRMSE;
5. factorized prediction beats the frozen-soma attacker in at least 8 / 12 holdouts;
6. no actual cluster run reaches the conservative -20 mV somatic spike guard.

If transport/oracle reconstruction passes but the frozen nonlinear operator
fails, classify:

~~~text
TRANSPORT_FACTORIZATION_VALID_LOCAL_OPERATOR_NOT_PORTABLE
~~~

If even the transport oracle fails:

~~~text
LINEAR_TRANSPORT_COMPOSITION_INADEQUATE
~~~

No kernel, stretch magnitude, branch, input dose or temporal alignment is tuned
after seeing the result.

## Interpretation

Passing would earn a genuinely reusable computational object:

~~~text
state/input-specific local nonlinear operator
                 x
geometry-specific transport operator
~~~

with geometry changes handled by replacing T rather than relearning N.

Failure is equally informative: it tells us exactly which attempted modular
boundary does not survive composition.

This is an architecture audit, not Gate 25.


## Implementation correction before scientific interpretation

The first executable attempt produced absurd original-geometry transport errors
of order 1e5 NRMSE. That cannot be a biological result because the operator
failed before any holdout.

The cause was identified before interpreting the architecture: the tiny 0.001
nA impulse response was baseline-corrected with a single pre-event voltage
number. Residual slow settling of the released model was therefore divided by
0.001 and mistaken for impulse response.

The corrected kernel is the causal matched-control difference:

~~~text
h(t) = [V_pulse(t) - V_no-pulse(t)] / 0.001 nA
~~~

and the event-driven soma trace is likewise:

~~~text
V_event(t) - V_no-event(t).
~~~

No scientific parameter or preregistered ruler changed:

- same six branches;
- same plus/minus 20 percent holdouts;
- same 24 virtual synapses;
- same 0.001 nA one-step impulse;
- same temporal alignment;
- same pass/fail thresholds.

The invalid run is retained as an implementation failure, not counted as a
scientific test of factorization.
