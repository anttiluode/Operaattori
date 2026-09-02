# Portable reduced circuit

This directory extracts the reusable mathematical object earned by the
Operaattori audits from the NEURON/FCI experiment harness.

`green_circuit_numpy.py` has **no NEURON dependency**. It implements:

~~~text
V_local = V_baseline + G * J

J = [g_AMPA + g_NMDA_raw * B_gamma(V_local)]
    * (E_rev - V_local)

V_soma = T * J
~~~

where:

- `G` is a measured local current-to-voltage Green/impulse matrix;
- `T` is a measured site-to-output transport kernel bank;
- the nonlinear law is explicit rather than a cached response waveform.

The module accepts any compatible measured kernels. It is not restricted to
cell 1125 or to three sites.

## Minimal shape contract

~~~python
from reduced.green_circuit_numpy import GreenCircuit

circuit = GreenCircuit(
    local_kernel_mV_per_nA_sample=G,   # [site, site, time]
    soma_kernel_mV_per_nA_sample=T,    # [site, time]
    baseline_local_mV=V0,              # [site, time]
)

result = circuit.solve(
    g_ampa_uS=gA,                      # [site, time]
    g_nmda_raw_uS=gN,                  # [site, time]
)
~~~

`result.soma_depolarization_mV` is the reduced output.

For delayed inputs, `shift_template` and `timed_conductance_matrix` provide
causal finite-window shifts without wraparound.

## Causal state-space runtime

The cross-cell nonlinear audit exposed a numerical boundary that the original
whole-waveform Green fixed point did not survive. `causal_graph_circuit.py`
extracts the passing organization as a simulator-free NumPy runtime:

~~~text
passive_(n+1) = P @ state_n
z             = passive_sites + R @ J(rest + z)
state_(n+1)   = passive_(n+1) + X @ J
soma          = state_(n+1)[soma_node]
~~~

The local implicit equation is solved with the analytic NMDA-current Jacobian
and damped Newton/backtracking. The module also contains a small dense reference
compiler from `G, C -> P, X` for tests and demos. The 24-cell scientific audit
used sparse LU instead; large morphologies should keep that sparse compilation
path rather than materializing a dense state matrix.

This distinction is deliberate: the earlier global waveform Picard arm
converged in only 207/288 cases, while the causal organization converged in all
288 archived comparisons.

~~~python
from reduced.causal_graph_circuit import CausalGraphCircuit

circuit = CausalGraphCircuit.from_dense_passive_graph(
    G_uS=G,
    C_nF=C,
    site_nodes=[i0, i1, i2],
    soma_node=soma,
    dt_ms=0.025,
)
result = circuit.run(g_ampa_uS=gA, g_nmda_raw_uS=gN)
~~~

See [../compiler.html](../compiler.html) for the browser reference specimen and
[../results/CAUSAL_NONLINEAR_GRAPH_DIAGNOSIS.md](../results/CAUSAL_NONLINEAR_GRAPH_DIAGNOSIS.md)
for the archived 288-case receipt.

## Exact operator / geometry tangents

The compiled causal object is now differentiable without an autodiff framework.
`operator_tangent.py` propagates first-order changes through the whole chain:

~~~text
(G, C)
  |
  v
(P, X)
  |
  v
implicit local NMDA
  |
  v
state trajectory
  |
  v
soma trace
~~~

For (A = G + D), (D = diag(C/dt)), (P=A^{-1}D), and (X=A^{-1}B),
the passive compiler tangent is analytic:

~~~text
dP = A^-1 (dD - dA P)
dX = -A^-1 dA X
~~~

At every nonlinear time step the site tangent reuses the implicit Newton
Jacobian:

~~~text
K  = I - R diag(dJ/dV)
dz = K^-1 [ d(passive_sites) + dR J ]
~~~

The implementation accepts a stack of geometry/control parameters, so a caller
can propagate several independent (dG/dtheta, dC/dtheta) directions in one
run. Conductance programs are held fixed.

The unit tests check the compiled (dP,dX) against centered finite
differences, then check the complete soma-trace tangent through the causal NMDA
loop against an independently recompiled finite difference. A two-parameter
path is checked as well.

This is a derivative of the **compiled circuit**. It does not assert that a
particular dendritic shape objective is biologically preferred. The next
scientific use is to connect these tangents to segment-level metric changes in
the real morphology compiler and use them for sensitivity/adversarial
validation.

## What is and is not shipped here

The runtime is generic NumPy code.

This directory intentionally does **not** redistribute the FCI cell-1125
mechanism bundle or a cell-specific kernel pack. The public Operaattori
receipts use a pinned released FCI model, but the reusable code and the model
data are kept separate while provenance/licensing of a redistributed derived
kernel pack is treated conservatively.

See:

- [../ARCHITECTURE.md](../ARCHITECTURE.md)
- [../results/GREEN_CIRCUIT_AUDIT.md](../results/GREEN_CIRCUIT_AUDIT.md)
- [../results/TEMPORAL_GREEN_CIRCUIT_AUDIT.md](../results/TEMPORAL_GREEN_CIRCUIT_AUDIT.md)
