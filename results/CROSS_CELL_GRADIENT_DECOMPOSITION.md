# Cross-cell transport / NMDA geometry-gradient panel

## Question

On human cell 1125, the exact morphology tangent separated cleanly into

~~~text
full geometry sensitivity
    =
frozen-current transport sensitivity
    +
voltage-dependent NMDA feedback sensitivity
~~~

and the nonlinear term could rotate or reverse the sign of a geometry
perturbation.

Does that architecture recur across different released morphologies, or was it
a property of one unusually convenient branch?

## Locked panel

- FCI commit: `75ad8b4d81a7f51bf888b30650c543592340db06`
- **24 released morphologies**
- one deterministic longest apical branch per morphology
- three nonlinear sites: x = 0.25 / 0.50 / 0.75
- two metric directions per site: local length / local diameter
- **6 geometry directions per cell per drive**
- drives: 0.5x / 1.0x / 2.0x
- **72 nonlinear operating points**
- **144 geometry directions per drive**
- matched passive morphology compiler; same causal local NMDA closure

Classification:

`CROSS_CELL_TRANSPORT_NMDA_REGIME_PANEL_VALID`

All 72 operating points converged.

## Main result

| drive | median max local depol | median full/transport cosine | median feedback/transport norm | cells feedback > transport | full-vs-transport sign reversals |
|---:|---:|---:|---:|---:|---:|
| 0.5x | 53.7 mV | 0.784 | 0.827 | 11 / 24 | 35 / 144 |
| 1.0x | 65.5 mV | 0.537 | 1.468 | 17 / 24 | 39 / 144 |
| 2.0x | 68.8 mV | 0.564 | 1.735 | 23 / 24 | 35 / 144 |

Across the complete panel:

~~~text
cells where feedback exceeds transport at >=1 drive      24 / 24
cells with >=1 full-vs-transport sign reversal           23 / 24
~~~

So the cell-1125 decomposition is not a single-cell curiosity.

As the local nonlinear sites are driven harder, the median feedback term grows
from smaller than transport to substantially larger than transport. The full
geometry-gradient vector correspondingly stops looking like the passive
transport gradient.

At 0.5x some cells are already in a strong nonlinear regime: the median maximum
local depolarization is 53.7 mV and 11/24 cells already have feedback norm
greater than transport norm. Therefore this panel does **not** establish a
universal "0.5x = passive regime" threshold. It establishes recurrence of the
decomposition across morphologies.

## The sign reversals are overwhelmingly diameter effects

Breakdown of full-vs-transport sign reversals:

| drive | length reversals | diameter reversals |
|---:|---:|---:|
| 0.5x | 4 | 31 |
| 1.0x | 1 | 38 |
| 2.0x | 0 | 35 |

The transport-only length derivative is negative at all 72 site positions at
1x and 2x, and the full system almost always preserves that sign.

Diameter is different. At 1x:

~~~text
transport-only positive diameter directions     31 / 72
full nonlinear positive diameter directions     51 / 72
diameter sign reversals                          38 / 72
~~~

At 2x:

~~~text
transport-only positive diameter directions     27 / 72
full nonlinear positive diameter directions     62 / 72
diameter sign reversals                          35 / 72
~~~

This makes the architecture distinction sharper:

> **Length sensitivity is comparatively transport-dominated in this site-level
> panel. Diameter sensitivity is strongly revalued by the local nonlinear
> state.**

That statement is about these three-site perturbations, not every possible
morphology direction.

## Heterogeneity matters

The aggregate trend is strong but individual morphologies differ.

At 0.5x the minimum full-vs-transport cosine is **-0.162**; at 1x it reaches
**-0.310**. In those cells, the nonlinear feedback does not merely amplify the
transport gradient — it rotates the six-dimensional geometry sensitivity far
enough that the two vectors point partly in opposite directions.

One morphology, rat L2/3 `229_5`, is the only cell with no sign reversal in
the three tested drive conditions. Its feedback/transport norm nevertheless
rises from 0.174 at 0.5x to 1.218 at 2x. So even the no-sign-flip case enters a
feedback-dominated magnitude regime.

A secondary species split is visible in this panel — rat cells have larger
median feedback/transport ratios than human cells at all three drives — but
that comparison was not the locked primary question and should not be treated
as a species claim without a dedicated audit.

## Architecture earned by the panel

The useful reusable statement is now:

~~~text
morphology metric theta
        |
        v
compiled transport operator T(theta)
        |
        v
local voltage state
        |
        v
voltage-dependent nonlinear current J(V)
        |
        +---- feedback into the same transport system
        |
        v
observable output
~~~

For a morphology perturbation,

~~~text
d output / d theta
    =
transport contribution
    +
local nonlinear feedback contribution
~~~

The transport term is determined by the compiled metric graph.

The nonlinear term is operating-point dependent and, across 23/24 released
morphologies, is strong enough in at least one tested condition to reverse the
sign of at least one local geometry direction.

This is more precise than saying "morphology is nonlinear" or "a neuron is a
nonlinear transfer function."

**Morphology supplies a transport operator. Local nonlinear state changes the
functional value of perturbing that operator.**

## Fence

This panel uses one deterministic apical branch per morphology and only the
three site-centered length/diameter directions.

It does not establish:

- a universal scalar threshold for the nonlinear regime;
- the same sign pattern on every branch;
- a biological growth objective;
- a learning rule;
- behavioral importance of the voltage changes;
- novelty over all prior differentiable morphology simulators.

Jaxley already demonstrates gradient-based optimization of compartment length
and radius. Operaattori's contribution here is the explicit causal
decomposition and its measurement across a zero-fit morphology-compiled panel.
