# Local Green matrix x synapse-law reduction

The transport-factorization audits passed, but they still supplied one measured
original-geometry nonlinear current waveform N_b(p) for every input pattern.

This audit removes that lookup.

## Reduced object

For each compact three-site branch and geometry g measure only the **linear**
impulse responses:

~~~text
G_g[i,j](t) : voltage at local site i from current injected at site j
T_g[j](t)   : soma voltage from current injected at site j
~~~

The nonlinear element is not fitted. It is the released synapse law itself:

~~~text
J_AMPA = g_AMPA(t) * (0 - V)

J_NMDA = g_NMDA_raw(t)
         * B_gamma(V)
         * (0 - V)

B_gamma(V) = 1 / (1 + exp(-gamma * V) / 3.57)
gamma = 0.078 / mV
~~~

Raw AMPA and voltage-independent NMDA conductance kinetics are measured once
from one original-geometry HUMAN synapse event, then reused everywhere. They
are mechanism kinetics, not a branch-response lookup.

The local voltages satisfy the nonlinear feedback equation:

~~~text
V_local(t) = V_noinput(t) + G_g * J(t)
J(t)       = synapse_law(V_local(t), input)
~~~

Solve this fixed point with locked damping 0.5 and relative current tolerance
1e-8. No held-out voltage/current/soma trace enters the solve.

Finally:

~~~text
Vhat_soma(t) = T_g * J_reduced(t)
~~~

## Panel

Reuse all six Gate-20 compact branches.

Geometries:

~~~text
1.00 x original  -- reconstruction control
0.80 x length    -- held out
1.20 x length    -- held out
~~~

Patterns:

~~~text
middle_single
outer_pair
triple
~~~

Multiplicity remains 8 per active site.

Total full-model comparison cases:

~~~text
6 branches x 3 geometries x 3 patterns = 54
~~~

## Attackers / diagnostics

For the 36 held-out geometry cases also compute:

1. frozen-soma attacker: original soma trace;
2. frozen-current factorization: T_g[J_original,p], the architecture that
   already passed the previous audit;
3. transport oracle: T_g[J_actual,g,p].

Thus the new reduced nonlinear circuit has to justify replacing the measured
original current lookup.

## Locked rulers

Earn:

~~~text
LOCAL_GREEN_MATRIX_X_SYNAPSE_LAW_REDUCES_RELEASED_NEURON
~~~

only if:

1. median transport-oracle soma NRMSE <= 0.02;
2. median reduced-circuit soma NRMSE across all 54 cases <= 0.05;
3. each pattern family's median reduced-circuit soma NRMSE <= 0.08;
4. median reduced current-waveform NRMSE <= 0.08;
5. on the 36 held-out geometry cases, reduced-circuit median soma NRMSE is
   <= 0.80 x the frozen-current-factorization median NRMSE;
6. reduced circuit beats frozen-current factorization in at least 24 / 36
   held-out cases;
7. every nonlinear fixed-point solve converges;
8. no actual comparison run crosses the -20 mV somatic spike guard.

If the transport oracle passes but the nonlinear circuit does not:

~~~text
TRANSPORT_REDUCTION_VALID_SYNAPSE_FEEDBACK_REDUCTION_INADEQUATE
~~~

No trace alignment, conductance gain, per-branch fit, damping change, pattern
selection or geometry selection is allowed after seeing results.

## Literature fence

This is not a claim to invent Green-function cable reduction, conductance-based
reduced neurons, two-layer dendritic subunits, or transfer-impedance-preserving
model reduction. Relevant prior art includes classical transient/Green
function cable methods and Neuron_Reduce (Amsalem et al., 2020,
doi:10.1038/s41467-019-13932-6).

The question here is only whether this exact released human-cell branch can be
causally decomposed into measured linear transport plus the released nonlinear
synapse law under held-out metric intervention.

This is an architecture audit, not Gate 25.


## Receipt — reduced circuit passes

The 54-case locked audit completed successfully.

~~~text
6 branches x 3 geometries x 3 patterns = 54

median transport-oracle soma NRMSE        0.0040
median reduced-circuit soma NRMSE         0.0043
median reduced current-waveform NRMSE     0.0038

held-out geometry only
  frozen-current factorization NRMSE      0.0282
  reduced-circuit soma NRMSE              0.0043
  reduced / frozen-current                0.1527
  reduced beats frozen-current            32 / 36

pattern-family reduced soma medians
  middle single                           0.0030
  outer pair                              0.0048
  triple                                  0.0070

fixed-point convergence                   54 / 54
actual soma spike guard                   0
~~~

Classification:

~~~text
LOCAL_GREEN_MATRIX_X_SYNAPSE_LAW_REDUCES_RELEASED_NEURON
~~~

This removes the per-pattern nonlinear-current lookup used by the previous
factorization tests.

The reduced object contains only:

1. the released AMPA conductance time course;
2. the released voltage-independent NMDA conductance time course;
3. the released HUMAN magnesium-block law with gamma = 0.078 / mV;
4. a measured 3 x 3 local current-to-voltage Green/impulse matrix for the three
   compact branch sites;
5. three measured site-to-soma transport kernels.

The site currents are solved self-consistently from the local voltages and the
released synapse law. No held-out soma, local-voltage, or nonlinear-current
trace is fitted.

The near equality between the transport oracle (0.40% median error) and the
fully reduced nonlinear circuit (0.43%) is the strongest part of the result.
At this operating regime, almost all of the full released model's response on
these assays is captured by **linear cable transport plus the explicit local
conductance nonlinearity**.

This is not a novelty claim for Green-function reduction or nonlinear dendritic
subunits. The result is the specific held-out causal reduction on the pinned
human cell-1125 reconstruction.

Compact CI receipt:
[results/operator_factorization/green_circuit_ci_summary.json](operator_factorization/green_circuit_ci_summary.json)

GitHub Actions: run 33494156490, job 99812336306.
