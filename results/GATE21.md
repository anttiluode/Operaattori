# Gate 21 — does the measured branch basis compute under a fixed input budget?

Gate 20 earned a specific abstraction in the released cell-1125 model:

    many semi-independent nonlinear branch subunits
                 -> one physical tree

Gate 21 asks whether that abstraction buys a compact input-output
factorization rather than merely a nicer description of the anatomy.

## Locked stimulus family

Use the exact same six long compact sections as Gates 16–20 and the same three
sites per section.

For every one of the 15 branch pairs, all six sites remain active. The global
budget is fixed at 48 virtual synapses, but the allocation between the two
branches changes:

    branch A   branch B   total
       12         36       48
       24         24       48
       36         12       48

Equivalently, the per-site multiplicities are 4+12, 8+8 and 12+4.

This keeps fixed, within a pair:

- the two branch identities;
- the same six physical synaptic sites;
- event time;
- total virtual synapse count;
- AMPA/NMDA mechanism and kinetics.

Only the allocation across the two branch compartments changes.

## Three rulers

For each redistribution, measure the actual somatic depolarization trace.

### 1. Total-input ruler

Because all three patterns have the same total input, this ruler predicts no
redistribution signature at all.

### 2. Independent-site ruler

Sum the six single-site soma traces at the exact site-specific dose.

This attacker knows which branches and sites were stimulated and preserves
their passive transfer and individual kinetics. It removes only interaction
between sites.

### 3. Nonlinear-subunit ruler

Measure each complete three-site branch alone at multiplicities 4, 8 and 12.
For a branch pair, predict the full-cell soma trace as:

    y_hat(A dose a, B dose b)
        = y_A_alone(a) + y_B_alone(b)

No branch-pair response is used by this ruler. It has only six branch transfer
curves and assumes cross-branch additivity.

## Primary metric: redistribution signature

For the three equal-budget patterns in each fixed branch pair, convert positive
soma AUC values to a centered log signature:

    s_k = log(AUC_k) - mean_j log(AUC_j)

This removes the pair's overall response scale and asks only whether the model
predicts **which redistribution is larger or smaller**.

Report RMSE of the signature for:

- zero signature (total-input ruler);
- independent-site prediction;
- nonlinear-subunit prediction.

Also report trace NRMSE for the two physical predictors and the actual
max/min AUC range produced by equal-budget redistribution.

## Preregistered positive classification

A HUMAN pass requires all of:

1. no conservative soma-spike contamination;
2. median equal-budget AUC range factor >= 1.05;
3. at least half of branch pairs have range factor >= 1.05;
4. median nonlinear-subunit signature error is at most half the
   independent-site error;
5. the nonlinear-subunit ruler beats the independent-site ruler on at least
   80% of branch pairs;
6. median nonlinear-subunit trace NRMSE <= 0.10.

The entire assay is repeated with Gate 20's sitewise rest-matched
gamma=0.062 attacker.

A pass earns only this statement:

> A small collection of measured nonlinear branch transfer functions
> factorizes equal-total-input somatic responses better than independent sites.

It does **not** earn learning, growth, optimality, or intelligence.

## Stopping line

If Gate 21 passes, Gate 22 must attack the exact dose lookup table: predict
held-out branch doses/redistributions from a smaller fitted/interpolated
subunit law. If Gate 21 fails, do not rescue it by scanning doses.


## Receipt — strong scalar signature, failed full-trace factorization

The locked assay ran on all 15 branch pairs with the corrected released-model
initialization.

~~~text
HUMAN

median equal-budget AUC range          1.3414x
pairs with range >= 1.05              15 / 15

signature RMSE
  total-input ruler                    0.1323
  independent-site ruler               0.1336
  nonlinear-branch ruler               0.0402

site / branch signature error          3.320x
branch ruler beats site               14 / 15
branch ruler beats total              15 / 15

full trace NRMSE
  independent-site ruler               0.7456
  nonlinear-branch ruler               0.1994

somatic spike guard                    0 / 15
~~~

The rest-matched gamma=0.062 control also produced redistribution sensitivity
and the branch ruler again predicted the centered AUC signature better than
independent sites:

~~~text
median equal-budget AUC range          1.2404x
site / branch signature error          2.764x
full branch-ruler trace NRMSE           0.3615
~~~

Classification:

~~~text
NO_CLEAN_BRANCH_SUBUNIT_FACTORISATION
~~~

This is deliberately a negative classification. The branch basis predicts
**which equal-budget redistribution produces more integrated somatic response**
surprisingly well, but it does not predict the complete somatic waveform within
the preregistered 10% error bound. HUMAN branch-ruler trace NRMSE is 19.94%.

So Gate 20's semi-independent compartments remain real, and Gate 21 shows that
input allocation across them matters, but the whole neuron is not cleanly the
sum of six measured nonlinear branch transfer functions.

The pre-registered stopping rule said not to scan doses to rescue a failure.
We therefore do not proceed to the held-out-dose Gate 22 that had been planned
for a positive factorization result.
