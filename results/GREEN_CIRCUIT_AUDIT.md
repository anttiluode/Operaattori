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
