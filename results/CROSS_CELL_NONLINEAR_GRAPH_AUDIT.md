# Cross-cell graph-generated nonlinear circuit audit

The passive cross-cell result now has a direct upstream map:

~~~text
morphology graph
      |
      v
hand-built passive cable solver
      |
      +--> local Green matrix G
      +--> soma transport T
~~~

Across 24 released morphologies that direct construction reached 0.21% median
joint G/T NRMSE without cross-cell fitting.

This audit asks whether those **graph-generated**, not NEURON-measured,
operators can drive the already-earned nonlinear AMPA/NMDA circuit on cells
outside the original cell-1125 work.

## Scientific question

Can the chain

~~~text
target morphology
      |
      v
graph-built G,T
      |
      v
explicit AMPA/NMDA voltage-feedback law
      |
      v
soma response
~~~

predict full nonlinear NEURON responses without first electrically measuring
the target cell?

## Target panel

Use all 24 pinned FCI morphologies under the same MATCHED_PASSIVE electrical
regime as the cross-cell passive audits.

For runtime while retaining whole-panel coverage, use the **three longest
apical sections** of every cell, selected deterministically by section cable
length with section-name tie break.

Each selected section uses three material addresses:

~~~text
x = 0.25, 0.50, 0.75
~~~

Cases:

~~~text
24 cells
x 3 branches
x 4 timing programs
= 288 full nonlinear comparisons
~~~

No cell is excluded, including human 2057 and rat L6 IPC.

## Matched passive cell

Exactly:

~~~text
Ra    = 150 ohm cm
Cm    = 1 uF / cm^2
Rm    = 20,000 ohm cm^2
E_pas = -70 mV

Na/Kv removed
~~~

This test is therefore about morphology-dependent cable loading plus a common
nonlinear probe, not species-specific fitted membrane physiology.

## Standardized nonlinear probe

Every material site receives the same released AMPANMDA_EMS mechanism with
explicitly locked parameters inherited from the successful cell-1125
Green-circuit work:

~~~text
AMPA event weight per active site  = 0.00088 uS x 8 = 0.00704 uS
NMDA / AMPA ratio                  = 0.00131 / 0.00088
gamma                              = 0.078 /mV
tau_r_AMPA                         = 0.20 ms
tau_d_AMPA                         = 1.70 ms
tau_r_NMDA                         = 0.29 ms
tau_d_NMDA                         = 43.0 ms
reversal                           = 0 mV
~~~

This is labelled HUMAN_PROBE. It is a controlled input law, not a claim that
the 24 released rat/human neurons share identical synaptic physiology.

The voltage-independent AMPA and raw-NMDA conductance templates are measured
**once** from the point-process mechanism under voltage clamp on a dummy
section. No target morphology contributes to those templates.

## Timing programs

Reuse the already-locked temporal family:

~~~text
synchronous   [ 0,  0,  0] ms
forward_5     [ 0,  5, 10] ms
reverse_5     [10,  5,  0] ms
spread_15     [ 0, 15, 30] ms
~~~

Reference event time:

~~~text
20 ms
~~~

Scored post-reference window:

~~~text
90 ms
~~~

Fixed time step:

~~~text
0.025 ms
~~~

## Reduced prediction

For each cell, build the full passive morphology graph exactly as in the
passing direct-cable audit.

At 0.025 ms, compute one morphology-derived impulse response bank for the nine
local source/target pairs and three soma transports of each selected branch.

No NEURON IClamp response is measured.

For one branch:

~~~text
v(t) = -70 mV + G_graph * J(t)

J_i(t) =
  [gA_i(t) + gN_raw_i(t) B_0.078(v_i(t))]
  [0 - v_i(t)]

y(t) = T_graph * J(t)
~~~

Solve the nonlinear current/voltage loop by the same damped fixed-point method:

~~~text
damping          0.5
relative tol     1e-8
max iterations   200
~~~

## Full-model reference

On the same matched-passive NEURON morphology, instantiate the same three
AMPANMDA_EMS point processes at x=.25/.50/.75 and deliver the same timing
program.

Record:

- three actual synaptic inward-current traces;
- soma voltage;
- local site voltages.

The full-model traces are used only for scoring.

## Diagnostics

### Graph transport oracle

Feed the **actual full-model synaptic currents** through graph-built T:

~~~text
y_oracle = T_graph[J_actual]
~~~

This diagnoses passive graph transport under the nonlinear drive without
granting the reduced model a measured target operator.

### Open-loop nonlinear attacker

Keep the same graph T and conductance templates but evaluate current at the
resting voltage only:

~~~text
J_open(t) = synapse_law(-70 mV, timed conductances)
y_open    = T_graph[J_open]
~~~

This tests whether the local voltage-feedback loop remains necessary
cross-cell.

## Locked primary ruler

Earn:

~~~text
MORPHOLOGY_GRAPH_X_NONLINEAR_LAW_PREDICTS_CROSS_CELL_RESPONSE
~~~

only if all hold:

1. median graph-transport-oracle soma NRMSE <= 0.02;
2. median reduced soma NRMSE <= 0.03;
3. median reduced current NRMSE <= 0.03;
4. every timing-family median reduced soma NRMSE <= 0.05;
5. median held-out-cell reduced soma NRMSE <= 0.04;
6. at least 20 / 24 cells have cell-median reduced soma NRMSE <= 0.10;
7. nonlinear fixed point converges in 288 / 288 cases;
8. reduced circuit beats the open-loop attacker in >= 90% of cases;
9. no target-cell impulse/impedance response enters G or T.

If the graph transport oracle passes but the nonlinear reduced current fails:

~~~text
CROSS_CELL_GRAPH_TRANSPORT_VALID_NONLINEAR_CLOSURE_NOT_PORTABLE
~~~

If graph transport itself fails:

~~~text
GRAPH_OPERATOR_NOT_ACCURATE_UNDER_CROSS_CELL_NONLINEAR_DRIVE
~~~

## Stronger descriptive label

If median reduced soma NRMSE <= 0.01 while the primary ruler passes, also
report:

~~~text
SUBPERCENT_CROSS_CELL_MORPHOLOGY_TO_NONLINEAR_RESPONSE
~~~

This is descriptive, not a separate novelty claim.

## Stopping rule

Do not tune:

- synaptic weight;
- gamma;
- timing delays;
- selected branches;
- passive constants;
- graph diameter correction;
- fixed-point damping;
- trace alignment;
- cell exclusions.

If the chain fails, diagnose whether error came from graph transport or local
nonlinear closure before any rescue.

## Prior-art fence

Passing would not make compartmental cable simulation, Green functions,
conductance synapses or reduced neuronal models novel.

The earned contribution would be the explicit held-out architecture test:
a target morphology is converted to G/T by an independent graph solver and
then to a nonlinear temporal response by the portable local feedback circuit,
without target-cell electrical calibration.

This is an architecture audit, not Gate 25.


## Receipt — global waveform fixed point fails portability

The locked 288-case audit completed.

~~~text
24 cells x 3 branches x 4 timing programs = 288

graph transport oracle soma NRMSE       0.0040
reduced soma NRMSE                      0.0042
reduced local-voltage NRMSE             0.0035
reduced current NRMSE                   0.0172

open-loop soma NRMSE                    0.5002
reduced / open-loop                     0.0083
reduced beats open-loop                 0.941

timing medians
  synchronous                           0.0041
  forward_5                             0.0044
  reverse_5                             0.0042
  spread_15                             0.0041

median cell soma NRMSE                  0.0049
cells with median soma NRMSE <= 0.10    21 / 24

rat median cell soma NRMSE              0.0134
human median cell soma NRMSE            0.0028

global fixed-point convergence          0.719
soma excursion guard                    0
~~~

Classification:

~~~text
CROSS_CELL_GRAPH_TRANSPORT_VALID_NONLINEAR_CLOSURE_NOT_PORTABLE
~~~

The classification is a **failure** because the preregistered global damped
waveform fixed point converged in only 207 / 288 cases.

### What survived

The graph-built passive transport survived nonlinear drive:

~~~text
actual full-model synaptic currents
        |
        v
graph-built T
        |
        v
soma

median error = 0.40%
~~~

Thus target-cell electrical operator measurement is not needed for downstream
transport on the median case.

The nonlinear feedback is also still essential. Removing local voltage feedback
while retaining the same conductance programs and graph transport gives 50.02%
median soma error.

### What failed

The failed cases are concentrated in a subset of morphologies/branches rather
than uniformly degrading the whole panel.

The worst cell is rat L6 UPC:

~~~text
median reduced soma NRMSE          2.4666
median current NRMSE              10.5507
median graph-transport oracle      0.0059
~~~

Rat L6 TPC also fails strongly:

~~~text
median reduced soma NRMSE          0.2508
median current NRMSE               0.5078
median graph-transport oracle      0.0032
~~~

Rat L6 IPC is a separate retained case where the hand-built passive graph itself
is less accurate:

~~~text
median reduced soma NRMSE          0.1119
median graph-transport oracle      0.0200
~~~

Across the 81 nonconverged fixed-point cases, the final relative-current update
has median approximately 0.28. Thirty-eight of those nonconverged cases still
happen to give <5% soma error, while others diverge catastrophically. Therefore
the convergence failure cannot be reclassified away by looking only at output
error.

### Diagnosis boundary

The current evidence separates:

~~~text
morphology graph -> passive transport        survives
global waveform Picard closure               not portable
~~~

It does **not** yet tell us whether the nonlinear physical closure itself fails
or whether solving an entire causal Volterra waveform as one damped Picard fixed
point is the wrong numerical algorithm for strongly loaded branches.

That is the next diagnosis. No damping change or probe retuning is opened.

GitHub Actions:
run 33537524016, job 99955378644.
