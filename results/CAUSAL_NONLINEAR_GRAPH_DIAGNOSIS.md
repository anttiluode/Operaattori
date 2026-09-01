# Causal nonlinear graph solver diagnosis

The preregistered cross-cell graph-generated nonlinear audit failed its primary
classification even though its median output error was very small:

~~~text
graph transport oracle            0.0040 soma NRMSE
global-waveform reduced circuit   0.0042
open loop                         0.5002

global Picard convergence         207 / 288
~~~

The graph transport therefore survived. The unresolved question is whether the
**nonlinear physical closure** is nonportable or whether solving an entire
causal Volterra waveform as one damped global Picard fixed point is the wrong
numerical algorithm for strongly loaded branches.

This diagnosis changes only the numerical solver.

## Same scientific panel

Reuse exactly:

~~~text
24 cells
3 longest apical branches per cell
3 sites x=.25,.50,.75
4 timing programs
288 cases

MATCHED_PASSIVE morphology
HUMAN_PROBE weight/gamma/kinetics
dt = 0.025 ms
90 ms scored window
~~~

No parameter from the failed run changes.

## Causal discrete graph equation

The morphology graph already supplies the passive compartment system:

~~~text
C dv/dt + G v = I
~~~

At one backward-Euler step:

~~~text
A = C/dt + G

v_(n+1)
  = A^-1 [ (C/dt) v_n + B J_(n+1) ]
~~~

where B injects the three synaptic currents at their graph compartments.

For a fixed branch, precompute:

~~~text
X = A^-1 B
R = S X
~~~

where S selects the three synaptic site voltages.

Let p be the passive prediction from the previous graph state:

~~~text
u = A^-1 (C/dt) v_n
p = S u
~~~

Then only the current time step is nonlinear:

~~~text
z = p + R J(-70 + z, gA_n, gN_n)
~~~

with z the three site depolarizations.

## Locked Newton solve

Solve the three-dimensional equation

~~~text
F(z) = z - p - R J(-70+z) = 0
~~~

using the analytic synapse derivative:

~~~text
B'(v) = gamma B(v) [1-B(v)]

dJ/dv =
  gN_raw B'(v) (0-v)
  - [gA + gN_raw B(v)]
~~~

Jacobian:

~~~text
dF/dz = I - R diag(dJ/dv)
~~~

Locked numerical settings:

~~~text
initial guess at first step    previous site depolarization (zero initially)
initial guess thereafter       previous solved site depolarization
Newton max iterations          30
absolute infinity-norm F tol   1e-10 mV
backtracking                   at most 8 halvings, accept first residual decrease
minimum line-search step       1/256
~~~

No damping parameter from the global Picard solver is reused or tuned.

After solving z:

~~~text
J = synapse_law(-70+z)
v_(n+1) = u + X J
~~~

and require the site entries of the full graph state to agree with z within
1e-8 mV.

## Simultaneous graph-transport oracle

For diagnosis, propagate a second passive graph state driven by the **actual
full-model synaptic currents**:

~~~text
v_oracle_(n+1)
 = A^-1 [ (C/dt) v_oracle_n + B J_actual_(n+1) ]
~~~

This is the causal state-space equivalent of the previous graph-T oracle.

## Full-model reference

Exactly the same full NEURON HUMAN_PROBE runs as the failed global-waveform
audit. No target-cell impulse or impedance measurement is introduced.

## Locked ruler

Earn:

~~~text
CAUSAL_MORPHOLOGY_GRAPH_NONLINEAR_CLOSURE_VALID
~~~

only if all hold:

1. median causal graph-oracle soma NRMSE <= 0.02;
2. median causal nonlinear soma NRMSE <= 0.03;
3. median causal current NRMSE <= 0.03;
4. median causal local-voltage NRMSE <= 0.03;
5. every timing-family median soma NRMSE <= 0.05;
6. median cell soma NRMSE <= 0.04;
7. at least 20 / 24 cells have cell-median soma NRMSE <= 0.10;
8. Newton converges at every time step of all 288 cases;
9. no soma excursion guard is crossed.

If this passes while the original global Picard convergence remains 207/288,
also classify the failure mechanism as:

~~~text
GLOBAL_WAVEFORM_PICARD_WAS_THE_NONPORTABLE_COMPONENT
~~~

If the causal graph oracle passes but the nonlinear solve still fails:

~~~text
CROSS_CELL_LOCAL_NONLINEAR_CLOSURE_ITSELF_FAILS
~~~

If the causal graph oracle fails:

~~~text
CROSS_CELL_PASSIVE_GRAPH_FAILS_UNDER_CAUSAL_DRIVE
~~~

If median nonlinear soma error <= 0.01 while the main ruler passes, report the
descriptive label:

~~~text
SUBPERCENT_CAUSAL_MORPHOLOGY_TO_NONLINEAR_RESPONSE
~~~

## Scope fence

Passing this diagnosis would support the full morphology-graph + nonlinear-law
architecture, but it would also show that the tiny FFT/Green-kernel runtime's
**global fixed-point algorithm** is not universally robust.

The causal solver keeps the full passive graph state. It is therefore a
scientific validation of the architecture, not yet a replacement for the
portable three-site NumPy runtime.

A later reduction would have to preserve this causal state-space stability.

## Stopping rule

Do not change the HUMAN_PROBE, graph discretization, dt, cell panel, branch
panel, delays, Newton tolerances, or line search after seeing the result.

This is a solver diagnosis, not Gate 25.


## Receipt — causal closure passes

The locked 288-case diagnosis passed.

~~~text
causal graph oracle soma NRMSE       0.00404
causal nonlinear soma NRMSE          0.00253
causal current NRMSE                 0.01722
causal local-voltage NRMSE           0.00080

timing medians
  synchronous                        0.00260
  forward_5                          0.00247
  reverse_5                          0.00251
  spread_15                          0.00242

median cell soma NRMSE               0.00214
cells with median soma NRMSE <= .10  23 / 24

Newton convergence                   288 / 288 cases
maximum Newton iterations            4
line-search failures                 0
site-consistency max                 < 1e-10 mV

rat median cell soma NRMSE           0.00724
human median cell soma NRMSE         0.00127
~~~

Classification:

~~~text
CAUSAL_MORPHOLOGY_GRAPH_NONLINEAR_CLOSURE_VALID
~~~

Failure mechanism of the previous audit:

~~~text
GLOBAL_WAVEFORM_PICARD_WAS_THE_NONPORTABLE_COMPONENT
~~~

Descriptive label:

~~~text
SUBPERCENT_CAUSAL_MORPHOLOGY_TO_NONLINEAR_RESPONSE
~~~

The scientific parameters are identical to the failed global-waveform run. The
difference is only the numerical organization of the nonlinear solve.

### The previous catastrophic cells recover

Rat L6 TPC:

~~~text
global-waveform cell median soma error      0.2508
causal cell median soma error               0.00024
~~~

Rat L6 UPC:

~~~text
global-waveform cell median soma error      2.4666
causal cell median soma error               0.00093
~~~

Those cells had excellent graph-transport oracles in the failed run, and the
causal solve now resolves them with three Newton iterations per step.

That is strong evidence that their earlier failures were solver divergence, not
a failure of the morphology-graph operator or the local HUMAN_PROBE law.

### Retained physical/discretization boundary

Rat L6 IPC remains the one cell outside the 10% cell-error fence:

~~~text
median causal soma NRMSE                    0.1119
median causal current NRMSE                 0.1572
median graph-oracle soma NRMSE              0.0200
all causal cases converged                  yes
~~~

This is the same morphology whose direct passive graph audit already had the
largest local-G error. The causal solver does not hide that earlier
discretization boundary.

### What is now supported

For the matched-passive 24-cell panel and the standardized HUMAN_PROBE:

~~~text
target morphology
      |
      v
hand-built passive cable graph
      |
      v
causal local voltage-feedback solve
      |
      v
nonlinear temporal soma response
~~~

requires no target-cell electrical calibration and reaches 0.25% median
soma-trace NRMSE across 288 cases.

The result does **not** establish species-specific synaptic physiology: the
same HUMAN_PROBE was intentionally applied to every morphology as a controlled
nonlinear test input.

### Runtime boundary

The earlier portable NumPy Green-kernel runtime solves the complete waveform by
global Picard iteration. That algorithm is no longer a universally supported
cross-cell runtime.

The robust cross-cell result uses causal state-space evolution of the passive
graph plus a three-dimensional implicit nonlinear solve at each time step.

A later reduced runtime must preserve that causal stability rather than merely
reusing the global waveform iteration.

Compact CI receipt:
[results/cross_cell_operator/causal_nonlinear_graph_ci_summary.json](cross_cell_operator/causal_nonlinear_graph_ci_summary.json)

GitHub Actions:
run 33538293022, job 99957910024.
