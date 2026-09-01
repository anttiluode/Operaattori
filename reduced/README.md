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
