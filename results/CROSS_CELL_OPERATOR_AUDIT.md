# Cross-cell morphology -> operator audit

The within-cell geometry interpolation result is deliberately easy in one
respect: every target belongs to the same reconstruction and lies between two
measured points on one scalar branch-length coordinate.

This audit removes that protection.

## Question

Can a neuron that was **never simulated during fitting** receive a useful local
Green/transport operator from morphology measurements alone?

The test uses all 24 released FCI morphologies at pinned commit:

~~~text
75ad8b4d81a7f51bf888b30650c543592340db06

12 rat
12 human
~~~

The authors' FCI values, species labels and synapse types are not predictor
features.

## Matched electrical regime

This is a morphology-to-operator test, not an exact reproduction of each
released cell's fitted somatic conductances.

Every morphology is instantiated with the released geometry and axon-stub
construction, then all Na/Kv mechanisms are removed and every section is reset
to the same passive parameters:

~~~text
Ra    = 150 ohm cm
Cm    = 1 uF / cm^2
Rm    = 20,000 ohm cm^2
E_pas = -70 mV
~~~

Thus the only cell-to-cell differences available to the electrical operator are
morphological.

This arm is labelled MATCHED_PASSIVE. It should not be confused with the
released FCI CODE_EXACT models, whose dendrites use Cm=2 and doubled leak while
soma/axon retain active Na/Kv.

## Deterministic branch panel

For each morphology:

1. collect apical sections;
2. sort by section cable length descending, tie-breaking by section name;
3. take the six longest sections;
4. define three material addresses at normalized x = 0.25, 0.50, 0.75.

No released synapse placement is used to select the branches.

This creates at most:

~~~text
24 cells x 6 branches = 144 operator examples
~~~

The held-out unit is always the **cell**, never an individual branch.

## Measured operator

At each branch, measure from one-dt +0.001 nA current impulses:

~~~text
G[i,j,t]   local voltage at site i from current at site j
T[j,t]     soma voltage from current at site j
~~~

A matched no-input trajectory is subtracted before dividing by the impulse
amplitude.

The scored window is locked to:

~~~text
dt          0.05 ms
impulse at  5 ms
post window 60 ms
~~~

No nonlinear synapse is involved in this first audit.

## Morphology-only features

No impedance, voltage, current, FCI, species or target-operator quantity enters
the predictor.

### Six local/tree features

For the selected apical section:

~~~text
1. log(section length)
2. log(mean diameter at x=.25,.5,.75)
3. log(path distance from soma to x=.5)
4. branch order from soma
5. log(total cable length of the section subtree)
6. log(1 + number of sections in the section subtree)
~~~

### Four whole-cell features

Computed directly from the instantiated morphology:

~~~text
7. log(total dendritic cable length)
8. log(total dendritic membrane area)
9. log(max dendritic path distance from soma)
10. log(total apical cable length)
~~~

These are available from an unseen morphology without measuring its electrical
operator.

## Predictor

Validation is leave-one-cell-out: 24 folds.

For every fold independently:

1. standardize features using the 23 training cells only;
2. build a training-only PCA basis separately for G and T;
3. keep the first 8 components (or fewer if rank-limited);
4. fit a fixed ridge map from standardized morphology features to PCA scores;
5. reconstruct the held-out cell's six G/T operator packs.

Ridge alpha is locked to 1.0. It is not cross-validated.

A second model uses only the six local/tree features. This tests whether global
cell morphology adds predictive information beyond the selected branch's own
cable context.

## Attackers and diagnostics

### Mean operator

Use the mean G/T pack of all training branches.

### Nearest morphology branch

In the same ten-dimensional standardized feature space, copy the operator from
the nearest training branch. The nearest branch can come from any training
cell, but never the held-out cell.

### PCA oracle

Project the actual held-out operator onto the **training-only** PCA basis and
reconstruct it.

This oracle is diagnostic only. It answers whether the cross-cell operator
family fits inside the locked low-dimensional representation before asking
morphology to predict coordinates within that representation.

## Metrics

For G and T separately:

~~~text
NRMSE = rms(predicted - actual) / rms(actual)
~~~

Report median branch error and median held-out-cell error.

Also report a joint error defined as the mean of the branch's G and T NRMSE.

## Locked primary ruler

Earn:

~~~text
CROSS_CELL_OPERATOR_PREDICTABLE_FROM_MORPHOLOGY
~~~

only if the ten-feature model satisfies all of:

1. median held-out branch joint NRMSE <= 0.20;
2. median held-out T NRMSE <= 0.20;
3. median held-out G NRMSE <= 0.20;
4. joint median error <= 0.80 x training-mean attacker;
5. joint median error <= 0.90 x nearest-morphology attacker;
6. morphology predictor beats nearest-morphology in >= 16 / 24 held-out cells;
7. no held-out cell is used in feature scaling, PCA or ridge fitting.

If the training-basis PCA oracle has joint median NRMSE <= 0.10 but morphology
prediction fails, classify:

~~~text
CROSS_CELL_OPERATOR_LOW_DIMENSIONAL_BUT_MORPHOLOGY_MAP_WEAK
~~~

If even the training-basis PCA oracle exceeds 0.10:

~~~text
CROSS_CELL_OPERATOR_FAMILY_NOT_CAPTURED_BY_LOCKED_BASIS
~~~

## Secondary questions

These do not affect the primary pass/fail decision:

- Does adding four global cell features improve over the six local/tree
  features?
- How do held-out human and rat errors compare?
- Does the same map extrapolate when trained only on rat cells and tested on
  human cells, and vice versa?

The species-transfer diagnostic uses the same fixed representation and ridge
settings. No retuning is permitted.

## Stopping rule

If this first morphology-only operator map fails, do **not** immediately rescue
it with neural networks, larger feature dictionaries, polynomial kernels or the
nonlinear NMDA circuit.

The failure mode must first be identified.

If it passes, the next justified assay is the expensive one: feed the predicted
held-out-cell G/T into the already-earned nonlinear Green circuit and compare it
against full NEURON responses on cells never used to fit the operator map.

## Prior-art fence

Green-function dendritic reduction, impedance matrices, reduced neuronal
models and nonlinear dendritic subunits are established prior art.

This audit does not claim otherwise.

It is also a direct descendant of the earlier
`Dig/AIZENBUD_GREENS_DIG.md` plan. Operaattori's contribution, if any, is the
held-out causal/operator portability experiment and the executable reduced
runtime, not invention of Green functions.

This is an architecture audit, not Gate 25.


## Receipt — cross-cell map fails

The full 24-cell leave-one-cell-out audit completed successfully.

~~~text
24 cells
6 deterministic apical sections per cell
144 measured branch operators

ten-feature morphology predictor
  median joint NRMSE                 0.3522
  median local-G NRMSE               0.2830
  median soma-T NRMSE                0.3736

attackers
  training-mean joint NRMSE          0.6686
  nearest-morphology joint NRMSE     0.3199

training-basis PCA oracle
  joint NRMSE                        0.0307

morphology / mean error              0.5268
morphology / nearest error           1.1010
held-out cells beating nearest       11 / 24

local/tree features only             0.3695
local + whole-cell features          0.3522
full / local-only                    0.9531
~~~

Classification:

~~~text
CROSS_CELL_OPERATOR_LOW_DIMENSIONAL_BUT_MORPHOLOGY_MAP_WEAK
~~~

The low PCA-oracle error is the central diagnostic. A held-out cell's operator
usually lies very close to the low-dimensional operator family spanned by the
other 23 cells, but the preregistered morphology coordinates do not locate it
accurately enough.

### Species diagnostics

~~~text
held-out rat cells, median cell joint NRMSE       0.2637
held-out human cells                              0.4263

train rat -> test human                           0.7842
train human -> test rat                           0.3804
~~~

The species labels themselves were never predictor features; these are only
post-run grouping diagnostics.

### Important retained failures

The most extreme cell is human L5 morphology 2057:

~~~text
morphology predictor       7.972
nearest branch             1.315
mean operator              9.021
PCA oracle                 0.091
~~~

It is also the largest morphology in this matched-passive panel by several
gross measures, including about 27,041 um total dendritic cable and a maximum
path over 2,100 um.

That makes it a genuine out-of-distribution stress case, not a reason to delete
the row.

For diagnosis only, removing 2057 after the fact does **not** rescue the
conclusion:

~~~text
23 remaining cells

morphology median joint NRMSE       0.3569
nearest median joint NRMSE          0.3366
PCA oracle                          0.0303
morphology beats nearest            11 / 23
~~~

So the failure is not one outlier.

### What the failure says

The operator representation itself survived the cross-cell jump much better
than the geometry-to-coordinate map:

~~~text
held-out operator -> training PCA family      ~3%
morphology features -> operator coordinates   ~35%
~~~

Adding the four whole-cell descriptors to the six local/tree descriptors only
improved median error from 36.95% to 35.22%.

The next justified question is therefore **not** "which larger regressor wins?"
It is whether the operator family needs a more physically appropriate
coordinate system than raw length/diameter/path/tree summaries.

No neural-network, polynomial-feature or nonlinear-NMDA rescue scan is opened.

GitHub Actions:
run 33506573377, job 99851921367.
