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
