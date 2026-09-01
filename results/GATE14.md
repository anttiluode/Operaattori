# Gate 14 — close a real branch through the released NMDA voltage gate

Gate 13 showed that passive bifurcations strongly load an otherwise identical
soma-to-tip path, but the passive transfer portfolio remained nearly rank one.

Gate 14 therefore adds the specific synaptic nonlinearity used by Aizenbud et
al., not an arbitrary activation function.

## Source-grounded mechanism

The paper uses conductance-based AMPA+NMDA excitation and the Jahr-Stevens
voltage-dependent magnesium block. The released FCI code gives the human
parameters used here:

    AMPA peak conductance      0.88 nS
    NMDA peak conductance      1.31 nS
    human gamma               0.078 /mV
    hybrid-B gamma            0.062 /mV
    [Mg2+]                    1 mM
    n                         1/3.57 mM^-1

Hybrid B deliberately keeps the human conductances and changes only gamma,
matching the paper's separation of synaptic strength from NMDA nonlinearity.

## Reduced assay

This gate is quasi-static at peak conductance. It does not pretend to reproduce
the full NEURON traces, synaptic rise/decay kinetics, soma spikes, or FCI.

For long real branches of cell 1125:

- derive the exact DC Green impedance matrix from the full passive branching tree;
- distribute 1, 2, 4, 8, 16, 32, and 48 simultaneous synapses along a branch;
- solve local voltages self-consistently because NMDA conductance depends on voltage;
- read the current arriving at a voltage-clamped soma boundary.

Attackers:

- fixed current: exact linear-superposition ruler;
- frozen NMDA block: same human AMPA/NMDA conductances but B(V) frozen at rest;
- hybrid B: same human conductances, smaller rat gamma;
- scalar geometry ruler: ask how much branch-to-branch enhancement is explained
  by median/max driving-point resistance plus branch length.

## What it can earn

A pass can establish branch-local nonlinear feedback on the matrix scaffold.
It cannot establish intelligence, learning, or reproduce the paper's FCI.

If the voltage-dependent condition materially exceeds the frozen-block ruler
and the effect is not just one scalar input-resistance variable, the next gate
earns the cost of a time-domain / NEURON-level synaptic assay.
