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
