# Electrotonic chart -> cross-cell operator audit

The first 24-cell cross-cell audit established:

~~~text
training-basis PCA oracle              0.0307 joint NRMSE
simple morphology -> operator          0.3522
nearest morphology branch              0.3199
~~~

So the held-out operator family is highly compressible, but raw
length/diameter/path/tree summaries are a poor coordinate chart.

This follow-up changes **only the predictor coordinates**.

## Physical hypothesis

Under the matched-passive regime, classical cable theory says the natural
length coordinate is electrotonic rather than metric:

~~~text
lambda(d) proportional to sqrt(d)
electrotonic length = integral ds / lambda(d)
characteristic admittance Y0 proportional to d^(3/2)
~~~

Branch-point loading is therefore naturally expressed through Rall-style
d^(3/2) combinations rather than section counts or raw cable length.

No voltage, current, impedance, FCI, species label, or target-operator
measurement is used as a predictor feature.

## Locked electrical constants

Exactly the previous matched-passive arm:

~~~text
Ra    = 150 ohm cm
Cm    = 1 uF / cm^2
Rm    = 20,000 ohm cm^2
E_pas = -70 mV
~~~

For diameter d:

~~~text
lambda_um(d) = 1e4 * sqrt((d_um * 1e-4) * Rm / (4 Ra))
~~~

## Same held-out object

Nothing else changes:

~~~text
24 released FCI morphologies
6 deterministic longest apical sections per cell
3 material sites x = 0.25, 0.50, 0.75
144 G/T operator packs

validation = leave one whole cell out
PCA components = 8
ridge alpha = 1.0
~~~

## Electrotonic feature chart

For each selected section, compute from morphology only:

1. section electrotonic length;
2. electrotonic soma-to-midpoint path length;
3. log characteristic admittance at midpoint, log(d_mid^(3/2));
4. taper, log(d_prox / d_dist);
5. distal Rall load ratio:
   sum child d_child^(3/2) / d_parent_dist^(3/2);
6. proximal Rall mismatch:
   d_parent^(3/2) /
   [d_selected_prox^(3/2) + sum sibling d_sibling^(3/2)];
7. total electrotonic cable length in the selected subtree;
8. maximum electrotonic depth from selected-section proximal end to any
   descendant endpoint;
9. subtree terminal admittance mass, sum_terminal d_terminal^(3/2);
10. subtree Rall-equivalent terminal diameter,
    [sum_terminal d_terminal^(3/2)]^(2/3);
11. total dendritic electrotonic cable length;
12. maximum soma-to-dendrite electrotonic path length.

For a terminal section the distal-load ratio is zero. For a soma-rooted section
the proximal mismatch uses neutral sentinel 1.0.

No species or layer code is included.

## Models compared

A. **Electrotonic chart** — the 12 features above. This is primary.

B. **Gross morphology chart** — the exact ten-feature predictor from the failed
audit.

C. **Combined chart** — all electrotonic + gross morphology coordinates.
Secondary only.

All use identical fold-wise standardization, training-only 8-D PCA and ridge
alpha 1.0.

## Attackers / diagnostics

Retain:

- training-mean operator;
- nearest training branch in the electrotonic feature space;
- training-only PCA oracle.

## Locked primary ruler

Earn:

~~~text
ELECTROTONIC_CHART_IMPROVES_CROSS_CELL_OPERATOR_MAP
~~~

only if the electrotonic model satisfies all:

1. median joint NRMSE <= 0.30;
2. median G NRMSE <= 0.30;
3. median T NRMSE <= 0.32;
4. joint median error <= 0.85 x previous gross-map error 0.3522;
5. joint median error <= 0.90 x electrotonic-nearest attacker;
6. beats electrotonic-nearest in >= 16 / 24 held-out cells;
7. PCA oracle remains <= 0.10 joint NRMSE.

If median joint NRMSE <= 0.20, classify more strongly:

~~~text
CROSS_CELL_OPERATOR_PREDICTABLE_FROM_ELECTROTONIC_MORPHOLOGY
~~~

If PCA remains good but the chart does not beat the gross map by the locked
amount:

~~~text
OPERATOR_LOW_DIMENSIONAL_ELECTROTONIC_CHART_STILL_INSUFFICIENT
~~~

## Secondary diagnostics

Report without altering pass/fail:

- combined electrotonic + gross chart;
- held-out rat vs human;
- rat -> human and human -> rat transfer;
- per-cell failures, retaining 2057;
- standardized ridge coefficient magnitudes into the operator PCA scores,
  averaged across folds.

No feature is selected or removed after seeing results.

## Stopping rule

Do not rescue failure with neural networks, polynomial features, target-cell
impedance measurements, species labels, deletion of 2057, direct PCA-coordinate
fits, or nonlinear-NMDA output scoring.

This is an architecture audit, not Gate 25.


## Receipt — electrotonic chart still insufficient

The locked 24-cell leave-one-cell-out audit completed.

~~~text
electrotonic chart
  median joint NRMSE                 0.4070
  median G NRMSE                     0.2643
  median T NRMSE                     0.4405

electrotonic nearest attacker        0.4068
training-basis PCA oracle            0.0307

electrotonic / previous gross        1.1556
electrotonic / nearest               1.0005
held-out cells beating nearest       14 / 24

gross morphology rerun               0.3522
combined gross + electrotonic        0.3062
~~~

Classification:

~~~text
OPERATOR_LOW_DIMENSIONAL_ELECTROTONIC_CHART_STILL_INSUFFICIENT
~~~

The primary physical chart therefore failed. Replacing metric morphology with
electrotonic/Rall summaries did not locate unseen cells more accurately.

### The useful secondary clue

The combined chart did improve over gross morphology alone:

~~~text
gross only joint NRMSE               0.3522
combined joint NRMSE                 0.3062
combined G NRMSE                     0.2323
combined T NRMSE                     0.3337
~~~

That is about a 13% relative reduction in median joint error, so the
electrotonic descriptors carry information that the gross morphology chart
misses. But they are complementary rather than sufficient.

The primary pass/fail decision remains a failure.

### Species diagnostics

~~~text
electrotonic LOCO rat median cell joint     0.3067
electrotonic LOCO human                     0.5577

rat -> human                                1.1211
human -> rat                                0.4211
~~~

The human side becomes harder, not easier, under the electrotonic-only chart.

### Largest coefficient diagnostics

Averaged across LOCO folds, the largest standardized ridge weights into the
operator PCA scores were:

~~~text
log midpoint characteristic admittance     15.52
section electrotonic length                  4.05
proximal/distal taper                        3.96
subtree total electrotonic length            3.85
soma-to-mid electrotonic path                3.67
distal Rall load ratio                       2.59
~~~

These weights are interpretive only; no feature was selected from them.

### Interpretation

The cross-cell operator family remains highly compressible:

~~~text
training-basis PCA oracle                    0.0307
~~~

but neither gross geometry nor a compact Rall/electrotonic summary is a good
coordinate chart.

The combined improvement suggests both carry real information, while the
remaining ~31% error points to missing **distributed boundary/loading
structure** rather than one more scalar feature.

The next justified question is therefore whether the morphology's full cable
graph itself is the appropriate coordinate: can a direct matched-passive graph
solver reconstruct G/T across the 24 cells without any cross-cell fitting?

No neural-net, polynomial-feature or target-impedance rescue is opened.

GitHub Actions:
run 33527487507, job 99921928756.
