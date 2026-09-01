# Gate 16a — can the released cell-1125 NEURON model run here?

Gate 15 exhausted the quasi-static reduction. Before building another reduced
mechanism, Gate 16a asks a software/representation question only:

> Can the authors' pinned cell-1125 model compile and execute its own released
> time-domain AMPA+NMDA mechanism inside this repository's CI?

This is deliberately **not** the scientific Gate 16 result.

## Pinned source

The CI checkout uses FCI commit:

    75ad8b4d81a7f51bf888b30650c543592340db06

and model:

    human/eyal/Human_L23_PC_0603_11_937_Eyal_passive_dends_simple_soma

with morphology:

    2013_03_06_cell11_1125_H41_06.asc

The model's own get_standard_model.py builds the passive-dendrite,
active-soma/axon cell and distributes the human synapses. Its own mods/
directory is compiled with nrnivmodl.

## Smoke event

After model construction, choose one deterministic dendritic synapse from the
released synapse table, schedule one manual event through its existing NetCon,
and record:

- somatic voltage;
- local dendritic voltage;
- AMPA and NMDA currents;
- AMPA and NMDA conductances.

The gate verifies the released human parameters (gamma=0.078, event weight
0.00088 uS) and requires finite, nonzero responses.

## What this earns

Only this:

    RELEASED_FCI_CELL1125_TIME_DOMAIN_RUNS

If it passes, Gate 16 may use the released kinetic model for a small matched
clustered-vs-dispersed temporal assay. It does not reproduce FCI and does not
unlock growth.
