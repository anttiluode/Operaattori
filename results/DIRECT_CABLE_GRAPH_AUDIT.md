# Direct morphology-graph passive operator audit

The cross-cell operator family is highly compressible, but two learned
morphology charts failed:

~~~text
gross morphology chart              0.3522 joint NRMSE
electrotonic chart                  0.4070
gross + electrotonic                0.3062
training-basis PCA oracle           0.0307
~~~

This audit stops searching for another low-dimensional input chart.

## Question

Can the **full loaded cable graph itself** generate the local Green matrix and
site-to-soma transport operator of an unseen morphology, with no cross-cell
fitting?

## Input

For each of the same 24 released FCI morphologies, instantiate exactly the same
MATCHED_PASSIVE geometry used by the two cross-cell audits.

The direct solver receives only:

- section topology;
- section length;
- segment diameter;
- segment membrane area;
- parent/child connection locations;
- fixed passive constants.

No NEURON voltage/current/impedance trace is used to construct the prediction.

NEURON is used only as the morphology importer and as the reference simulator.

## Locked passive constants

~~~text
Ra    = 150 ohm cm
Cm    = 1 uF / cm^2
Rm    = 20,000 ohm cm^2
E_pas = -70 mV
~~~

All active Na/Kv mechanisms are removed, as in the previous cross-cell audits.

## Discrete cable graph

Each NEURON segment center becomes one membrane compartment.

For compartment i:

~~~text
C_i = Cm * membrane_area_i
gL_i = area_i / Rm
~~~

with units converted to nF and uS.

Each half-segment is treated as a cylindrical axial resistor:

~~~text
R_half
  = Ra * (L_segment/2) / [pi (d/2)^2]
~~~

using segment-center diameter.

Segment boundaries and branch junctions are represented as zero-capacitance
nodes. They are eliminated analytically.

For a junction with half-segment conductances g_k connected to membrane
compartments k, elimination contributes the Schur-complement Laplacian:

~~~text
L_junction
  = diag(g) - g g^T / sum(g)
~~~

This correctly couples siblings through a common branch point rather than
pretending every child has an independent parent edge.

The complete matched-passive system is:

~~~text
C dv/dt + (G_leak + L_axial) v = I(t)
~~~

for depolarization relative to -70 mV.

## Time stepping

Use backward Euler at the exact reference step:

~~~text
dt = 0.05 ms
~~~

For each step:

~~~text
[C/dt + G] v_(n+1)
  = [C/dt] v_n + I_(n+1)
~~~

One sparse LU factorization is reused for all sources and all time steps of one
cell.

No trace alignment or fitted temporal correction is permitted.

## Same scored operator

Exactly the previous panel:

~~~text
24 cells
6 deterministic longest apical sections per cell
3 sites at x = 0.25, 0.50, 0.75
144 branch operator packs

impulse amplitude     +0.001 nA
impulse duration      0.05 ms
post window           60 ms
~~~

For every branch score:

~~~text
G[3,3,t]   local current -> local voltage
T[3,t]     current -> soma voltage
~~~

against the NEURON MATCHED_PASSIVE reference.

## No learning

There is:

- no training/test split because there are no fitted cross-cell parameters;
- no PCA;
- no ridge regression;
- no species label;
- no FCI value;
- no target-cell calibration;
- no target impedance measurement.

Every cell is predicted independently from its own morphology graph plus the
same fixed physical constants.

## Locked rulers

Earn:

~~~text
MORPHOLOGY_GRAPH_GENERATES_PASSIVE_OPERATOR
~~~

only if all hold:

1. median joint G/T NRMSE <= 0.05;
2. median G NRMSE <= 0.05;
3. median T NRMSE <= 0.05;
4. median held-out-cell joint NRMSE <= 0.06;
5. at least 20 / 24 cells have median joint NRMSE <= 0.10;
6. worst non-pathological cell joint NRMSE is reported, not hidden;
7. no trace alignment, fitted gain, cell-specific correction, or segment
   selection is introduced after seeing the result.

A looser but still useful classification is earned if median joint NRMSE <=
0.10:

~~~text
MORPHOLOGY_GRAPH_APPROXIMATES_PASSIVE_OPERATOR
~~~

If median error exceeds 0.10:

~~~text
HAND_BUILT_CABLE_DISCRETIZATION_INADEQUATE
~~~

## Interpretation fence

Passing would **not** be a novelty claim for cable theory, compartmental
models, Green functions, or passive neuronal simulation. Those are classical.

Its value here would be architectural:

~~~text
morphology graph
   |
   v
direct physical operator construction
   |
   v
G, T
   |
   v
already-earned nonlinear reduced circuit
~~~

That would explain why scalar cross-cell charts failed even though the operator
family itself was low-dimensional: the correct upstream representation retains
the loaded tree until after the cable physics has acted.

## Stopping rule

If the graph solver fails, diagnose discretization error before altering the
biology. Do not rescue it with learned corrections.

If it passes, the next test may feed graph-generated G/T into the portable
nonlinear Green circuit on cells never used in the original cell-1125 work.

This is an architecture audit, not Gate 25.


## Receipt — morphology graph generates the passive operator

The locked 24-cell direct-physics audit passed.

~~~text
24 cells
6 deterministic apical sections per cell
144 branch operator packs

median joint G/T NRMSE             0.0021
median local G NRMSE               0.0013
median soma T NRMSE                0.0024

median held-out-cell joint NRMSE   0.0018
cells with cell median <= 0.10     23 / 24
~~~

Classification:

~~~text
MORPHOLOGY_GRAPH_GENERATES_PASSIVE_OPERATOR
~~~

No cross-cell fitting, PCA coordinates, species labels, impedance measurements
or target-cell calibration were used.

The direct solver received the instantiated morphology graph, fixed passive
constants and the same material addresses, then assembled membrane
capacitances, leaks and axial conductances. Zero-capacitance branch junctions
were eliminated with the preregistered Schur-complement construction.

### The former "impossible" cell

Human L5 morphology 2057 was the extreme failure of the learned morphology
chart:

~~~text
gross morphology predictor       about 7.97 joint NRMSE
nearest training branch          about 1.32
training PCA oracle              about 0.09
~~~

The direct morphology graph gives a cell-level median joint error of roughly:

~~~text
0.0001
~~~

So 2057 was not intrinsically difficult for passive operator construction. It
was difficult to compress into a small feature vector.

### Retained failure mode

The worst cell is rat L6 IPC:

~~~text
median G NRMSE                   0.1761
median T NRMSE                   0.0348
median joint NRMSE               0.1101
~~~

It is retained.

The error is concentrated much more strongly in local G than soma T, consistent
with the hand-built cylindrical half-segment approximation being less faithful
to some local diameter/3-D discretizations than to downstream somatic
transport.

No post-hoc diameter correction, NEURON axial-resistance query or fitted gain
was introduced.

### Architectural interpretation

The failed learned charts and the successful direct solver now separate three
things:

~~~text
held-out cross-cell operator family
    is low-dimensional                 yes (~3% PCA oracle)

small scalar morphology chart
    locates an unseen cell             no  (~35-41%)

full loaded cable graph
    generates the passive operator     yes (~0.2%)
~~~

The conclusion is not that cable theory is new. The conclusion is that
Operaattori's upstream representation should preserve the loaded morphology
graph until after the physical cable operator is constructed.

The current cross-cell architecture is therefore:

~~~text
morphology graph
      |
      v
passive cable construction
      |
      +--> local Green matrix G
      |
      +--> site-to-soma transport T
~~~

The nonlinear AMPA/NMDA circuit has not yet been re-run on this cross-cell
graph-generated operator panel. That remains a separate scientific test.

Compact CI receipt:
[results/cross_cell_operator/direct_cable_graph_ci_summary.json](cross_cell_operator/direct_cable_graph_ci_summary.json)

GitHub Actions:
run 33528552446, job 99925516858.
