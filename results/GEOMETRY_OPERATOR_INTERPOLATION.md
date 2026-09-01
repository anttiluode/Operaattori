# Geometry -> operator interpolation audit

The reduced Green circuit now generalizes across branch identity, intrinsic
metric perturbation, input pattern and event timing. One important dependency
remains: every geometry-specific Green/transport operator has so far been
measured directly.

This audit asks whether the operator itself lies on a simple geometry-indexed
family.

## Endpoint operators

For every one of the six Gate-20 compact branches, measure the complete reduced
linear pack only at:

~~~text
lambda = 0.80
lambda = 1.20
~~~

where lambda multiplies only the selected branch section's intrinsic cable
length.

Each endpoint pack contains:

~~~text
V0_lambda(t)       no-input local voltage at the 3 sites
G_lambda[i,j](t)   3x3 local current -> voltage kernels
T_lambda[j](t)     3 site -> soma transport kernels
~~~

The HUMAN conductance templates and nonlinear magnesium-block law are the same
ones already earned by the Green-circuit audits.

## Held-out target geometries

Score, but do not use to construct the predicted operator:

~~~text
lambda* = 0.90
lambda* = 1.00
lambda* = 1.10
~~~

For a target lambda*, define

~~~text
alpha = (lambda* - 0.80) / (1.20 - 0.80)

V0_hat = (1-alpha) V0_0.8 + alpha V0_1.2
G_hat  = (1-alpha) G_0.8  + alpha G_1.2
T_hat  = (1-alpha) T_0.8  + alpha T_1.2
~~~

No coefficient is fitted from target traces.

## Temporal panel

Reuse the locked temporal programs:

~~~text
synchronous   [ 0,  0,  0] ms
forward_5     [ 0,  5, 10] ms
reverse_5     [10,  5,  0] ms
spread_15     [ 0, 15, 30] ms
~~~

The predicted target response is obtained only by solving the released
AMPA/NMDA law inside the interpolated operator pack.

## Cases

~~~text
6 branches
x 3 held-out target geometries
x 4 temporal programs
= 72 held-out target responses
~~~

## Attacker

Use the nearest endpoint pack without interpolation:

~~~text
lambda*=0.90 -> use 0.80 pack
lambda*=1.00 -> use 0.80 pack (locked lower-end tie break)
lambda*=1.10 -> use 1.20 pack
~~~

The same synapse law and timed conductances are used. This asks whether
interpolating the operator adds predictive value beyond merely choosing a
nearby measured geometry.

## Oracles

For diagnosis only, target geometries also receive directly measured target
G/T kernels after the prediction is fixed.

Report:

1. target transport oracle T_target[J_actual];
2. target reduced-circuit oracle using target V0/G/T plus the same synapse law;
3. interpolated-operator reduced prediction;
4. nearest-endpoint reduced attacker.

## Locked rulers

Earn:

~~~text
GEOMETRY_INTERPOLATES_REDUCED_OPERATOR
~~~

only if:

1. median target transport-oracle soma NRMSE <= 0.01;
2. median target reduced-oracle soma NRMSE <= 0.02;
3. median interpolated-operator soma NRMSE <= 0.02;
4. each target-scale median interpolated soma NRMSE <= 0.03;
5. each timing-family median interpolated soma NRMSE <= 0.03;
6. median interpolated current-waveform NRMSE <= 0.02;
7. interpolated median soma NRMSE <= 0.75 x nearest-endpoint median NRMSE;
8. interpolation beats the nearest endpoint in at least 48 / 72 cases;
9. all interpolated nonlinear fixed-point solves converge;
10. no actual target run crosses the -20 mV soma guard.

If target direct operators work but interpolation fails:

~~~text
REDUCED_CIRCUIT_VALID_OPERATOR_GEOMETRY_MAP_NOT_LINEAR
~~~

If target direct operators themselves fail:

~~~text
TARGET_GEOMETRY_REDUCTION_FAILED
~~~

No polynomial fit, alignment, branch weighting, endpoint movement, target
selection, timing change or conductance correction is allowed after seeing the
result.

This is an architecture audit, not Gate 25.
