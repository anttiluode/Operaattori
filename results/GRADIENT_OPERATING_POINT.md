# Geometry gradient operating-point audit — cell 1125

## Why this audit exists

The real cell-1125 metric tangent is valid: analytic derivatives through the
direct cable compiler and causal NMDA closure match centered metric
recompilation.

That did **not** answer two harder questions:

1. Is the spatial sign pattern stable, or was the x≈0.25 length sign flip a
   single operating-point accident near the NMDA knee?
2. Is the gradient large and linear enough to move the output under finite
   geometry changes, or is it mathematically correct but practically inert?

This audit keeps the same human cell, same apical branch, same three nonlinear
sites and same 0/5/10 ms timing program, then changes only the drive scale.

## Locked object

- FCI commit: `75ad8b4d81a7f51bf888b30650c543592340db06`
- human L2/3 morphology: `1125`
- full compiled graph: **1653 compartments**
- tested branch: `PassiveDendsSimpleSomaModel[0].apic[100]`
- branch compartments: **31**
- metric directions per operating point: **62**
- sites: x = 0.25 / 0.50 / 0.75
- drive scales: 0.25 / 0.50 / 0.75 / 1.00 / 1.50 / 2.00 / 3.00
- finite responses: fresh full `G,C` recompilation, not tangent-only estimates

Classification:

`GRADIENT_OPERATING_POINT_AUDIT_VALID`

## Result 1 — the pretty baseline sign flip is not universal

| drive | soma peak mV | max local depol mV | length + / − | diameter + / − | first negative length x | gradient L2 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.4597 | 12.95 | 0 / 31 | 10 / 21 | 0.016 | 0.0369 |
| 0.50 | 0.9198 | 24.29 | 0 / 31 | 5 / 26 | 0.016 | 0.0712 |
| 0.75 | 1.4798 | 35.30 | 0 / 31 | 0 / 31 | 0.016 | 0.1424 |
| **1.00** | **2.3557** | **46.65** | **7 / 24** | **0 / 31** | **0.242** | **0.3501** |
| 1.50 | 3.9380 | 60.57 | 0 / 31 | 19 / 12 | 0.016 | 0.3092 |
| 2.00 | 4.6616 | 64.90 | 0 / 31 | 23 / 8 | 0.016 | 0.5256 |
| 3.00 | 5.3226 | 67.59 | 0 / 31 | 24 / 7 | 0.016 | 0.7262 |

At the released 1× program, local length is slightly positive in the seven
compartments before the first nonlinear site and negative from x≈0.242 onward.

That structure disappears on **both sides** of the operating point:

- at 0.25×, 0.50× and 0.75×, length is negative in 31/31 compartments;
- at 1.50×, 2× and 3×, length is again negative in 31/31 compartments;
- diameter is 31/31 negative at 0.75× and 1×, but flips to majority-positive at
  high drive (19/31, 23/31, 24/31 positive).

Therefore:

> **The baseline x≈0.25 sign flip is a nonlinear operating-point signature, not
> a general growth law.**

That is the main result of this audit.

## Result 2 — the analytic gradient remains useful at finite geometry budget

At drive = 1×, take the full 62-dimensional gradient of soma peak and move in
its steepest-ascent or steepest-descent direction under an L2 log-metric
budget. Then rebuild the full graph and rerun the causal nonlinear circuit.

### Steepest ascent

| L2 log budget | largest local change | predicted fixed-time Δ | actual peak Δ | peak change |
|---:|---:|---:|---:|---:|
| 0.01 | ~0.279% | +0.003501 mV | +0.003504 mV | +0.149% |
| 0.05 | ~1.40% | +0.017505 mV | +0.017533 mV | +0.744% |
| 0.10 | ~2.82% | +0.035009 mV | +0.035072 mV | +1.489% |
| 0.20 | ~5.72% | +0.070019 mV | +0.069994 mV | +2.971% |

At the largest tested budget, the fixed-time nonlinear response remains within
about **2.3%** of the first-order prediction.

The descent arm is similarly well behaved: −2.94% actual peak change at the
same 0.20 L2 budget.

So the gradient is not merely a correct infinitesimal. It remains quantitatively
useful over several-percent local metric changes on this branch.

This is a mathematical optimizer statement. The ascent direction is free to
shorten some cable and narrow some compartments. It is **not** a biological
growth rule.

## Result 3 — literal local lengthening is much weaker

To keep the growth-like question separate, two single-compartment
**lengthening-only** sweeps were run at drive = 1×.

### Before the first nonlinear site: strongest positive length gradient

Position: x = 0.016.

| length increase | actual soma-peak change | percent of peak |
|---:|---:|---:|
| +1.01% | +0.0000158 mV | +0.00067% |
| +5.13% | +0.0000796 mV | +0.00338% |
| +10.52% | +0.0001606 mV | +0.00682% |
| +22.14% | +0.0003221 mV | +0.01367% |

The warm pre-site growth effect is real but tiny.

### Downstream: strongest negative length gradient

Position: x = 0.661.

| length increase | actual soma-peak change | percent of peak |
|---:|---:|---:|
| +1.01% | −0.0004596 mV | −0.0195% |
| +5.13% | −0.0023447 mV | −0.0995% |
| +10.52% | −0.0048097 mV | −0.204% |
| +22.14% | −0.0101246 mV | −0.430% |

The downstream cooling effect is larger, but still modest.

So the safe interpretation is:

> **A useful distributed metric optimization signal exists. Literal growth-only
> effects on one segment are much smaller.**

That is not a failure. It prevents the differentiability result from being
inflated into a strong growth-computation claim.

## Prior-art correction

The repo must not claim that differentiable dendritic shape is unprecedented.

Jaxley (Nature Methods, 2025) explicitly supports gradients with respect to
morphological parameters and demonstrates gradient-based training of
compartment **length and radius** in a nonlinear single-neuron task.

Operaattori's result is narrower:

- explicit morphology → passive cable construction with zero fitting;
- hand-derived analytic metric tangent through junction elimination;
- local implicit NMDA tangent rather than generic autodiff through the whole
  simulator;
- validation on a released human reconstruction;
- exact-zero pose direction by construction.

Whether that narrower package is novel requires a separate literature and code
survey.

## Fence

This audit establishes operating-point dependence and finite-budget behavior
for one fixed-topology branch and one three-site timing family.

It does **not** establish:

- a biological growth objective;
- gradient descent by real dendrites;
- topology-changing growth;
- a behaviorally meaningful voltage threshold;
- cross-cell stability of the gradient field;
- novelty versus all prior differentiable neuron simulators.

The next scientific question is no longer “is the gradient real?” It is:

> **Does the operating-point-dependent geometry field generalize across branches
> and cells, and can one predict its regime changes from local cable/NMDA state
> rather than recomputing the full gradient map?**
