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


## Receipt — electrical compartment yes, NMDA locality advantage no

Gate 15 was run twice because the first definition of "clustered" was too
generous.

### Gate 15a — same unbranched run, spread over the whole branch

Eight long maximal unbranched runs were compared with passively matched sites
distributed over six to eight other branch runs.

~~~text
median passive Rinput match factor         1.069x
median passive soma-transfer match         1.039x
minimum dispersed branch runs              6

median clustered/dispersed Green coupling  4.915x

high-dose locality index:
  soma current (human/frozen)              0.9903
  local voltage (human/frozen)             0.9403
  max local voltage (human/frozen)         0.9620

branches with soma-current locality >1.05  0 / 8
branches with local-V locality >1.05       0 / 8
~~~

The sites called "clustered" here spanned roughly 293–385 um. That is one
unbranched dendritic run, but not necessarily one compact integration
compartment, so this negative result was not taken as the stopping result.

### Gate 15b — compact midpoint window

The same protocol was repeated with eight clustered reconstruction nodes
restricted to a 50-um physical window. Actual occupied spans were 39.4–48.5 um.

~~~text
median passive Rinput match factor         1.059x
median passive soma-transfer match         1.080x
minimum dispersed branch runs              6

median clustered/dispersed Green coupling  9.557x

high-dose locality index:
  soma current (human/frozen)              0.9241
  local voltage (human/frozen)             0.9196
  max local voltage (human/frozen)         0.9652

median local voltage at 48 synapses:
  compact HUMAN                            -2.61 mV
  compact FROZEN                           -5.43 mV
  dispersed HUMAN                         -11.14 mV
  dispersed FROZEN                        -17.93 mV

branches with soma-current locality >1.05  0 / 8
branches with local-V locality >1.05       0 / 8

HUMAN / hybrid-B locality                  1.0078
maximum nonlinear-solver residual          9.698e-10 mV
~~~

Classification:

~~~text
NO_ROBUST_NMDA_LOCALITY_ADVANTAGE
~~~

The passive electrical compartment is unambiguous: compact same-branch sites
are almost ten times more mutually coupled than the matched dispersed sites.
What fails is the stronger hypothesis that this locality specifically magnifies
the *human NMDA voltage-dependence* in the quasi-static assay.

The compact group is already pushed very near the 0-mV excitatory reversal
potential even with the magnesium block frozen. Human NMDA depolarizes it
further, but the relative human/frozen advantage is smaller there than for the
dispersed matched inputs. Measuring local voltage instead of soma clamp current
does not rescue the claim.

This does **not** say that dendrites lack compartments. It says that the
particular reduced mechanism tested here — compact branch locality selectively
amplifying the human NMDA voltage gate — is not supported.

## Revised stopping line

Do not rescue Gate 15 by adding growth or another invented nonlinearity.

If this line continues, the next justified escalation is a small time-domain
assay in the authors' released NEURON model for cell 1125, using its actual
AMPA/NMDA rise/decay kinetics and the same matched clustered/dispersed logic.
That would test whether the equilibrium approximation erased a temporal
compartment effect. If the released model is not practical to run in a bounded
gate, stop this mechanistic branch rather than replacing it with a homemade
surrogate.
