# Gate 15 — does one branch become a nonlinear compartment?

Gate 14 showed a modest effect of the released human NMDA voltage dependence,
but essentially no new dose-response dimension.

Gate 15 asks whether the missing ingredient is **where coactive synapses share
the dendritic tree**.

Aizenbud et al. argue that larger, more elaborate morphology likely increases
dendritic compartmentalization and semi-independent integration, while stronger
voltage-dependent NMDA supports nonlinear interactions among coactive
excitatory synapses. This gate is a reduced mechanistic probe of that proposed
bridge; it is not an experiment reported in their paper and is not an FCI
reproduction.

## Matched contrast

For each long maximal unbranched run in cell 1125:

1. choose eight distinct reconstruction nodes distributed across that one run;
2. for every clustered node, find a unique site on another branch with similar
   **driving-point impedance** and **absolute soma current transfer**;
3. force the dispersed control to occupy at least six distinct branch runs;
4. activate the same number of human AMPA+NMDA synapses in each arrangement.

Each selected site receives 1, 2, 4, then 6 identical synapses, so the eight-site
assay reaches 8, 16, 32, and 48 simultaneous synapses without changing the
matched site identities.

## Exact passive ruler

The all-node driving-point impedance is computed by exact O(N) rerooting of the
passive tree. Pairwise Green impedance matrices are then computed separately
for the clustered and dispersed sites.

Matching is done **before** seeing any nonlinear response, using only:

    log(Rinput)
    log(|Tsoma|)

A poor match is allowed to classify the gate as inconclusive.

## Nonlinear contrast

For each arrangement compute:

    R = soma response with voltage-dependent human NMDA
        ------------------------------------------------
        soma response with the same AMPA/NMDA conductances
        but magnesium block frozen at rest

The compartment-locality index is:

    L = Rclustered / Rdispersed

If L is robustly above one after passive matching, co-locating matched inputs on
one branch specifically amplifies the voltage-dependent NMDA contribution.

Hybrid B is also included: it keeps the human conductances but lowers gamma to
the rat value, providing a second ruler for the steepness of the voltage gate.

## Why this is harder than Gate 14

Gate 14 could be explained partly by a branch having a high local input
resistance. Gate 15 attacks that explanation directly by matching each
clustered site to a dispersed site with similar individual passive electrical
properties.

The arrangements should nevertheless differ in their **off-diagonal Green
coupling**: same-branch sites should perturb one another more strongly than
passively matched sites spread across the tree. That interaction matrix is the
candidate physical compartment.

## Outcome classes

    PASSIVE_MATCH_INADEQUATE
    NO_ROBUST_NMDA_LOCALITY_ADVANTAGE
    NMDA_LOCALITY_ADVANTAGE_PRESENT

Scientific failure is allowed to pass CI. Numerical failures are not.

## Stopping line

A positive Gate 15 result would earn a time-domain synaptic assay with the
paper's rise/decay kinetics. It would still not earn growth, intelligence, or an
FCI claim.

A negative result means the reduced quasi-static scaffold has hit another
ceiling and should not be rescued by adding developmental plasticity.
