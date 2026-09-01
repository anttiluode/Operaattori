# Gate 20 — are branches actually semi-independent nonlinear subunits?

Gate 19 killed the simple sequence-order story. Gates 16–18 still leave a
cleaner possibility: compact dendritic branches may behave as nonlinear
**compartments** rather than sequence processors.

Gate 20 tests the word `compartment` directly.

## Hierarchical superposition

Use the same six long compact cell-1125 sections selected before any nonlinear
outcome in Gates 16–19.

At a locked 8 virtual synapses per site, one active branch contains:

    3 sites x 8 = 24 synapses

and a branch pair contains 48 synapses total, keeping the global drive below
the high-dose pair that could invite somatic spikes.

### Within one branch

For each branch and condition:

    I_within = response(3 sites together)
               /
               sum(response(each site alone))

using local-voltage AUC at that branch's three sites.

### Between two branches

For all 15 pairs of the six branches:

    I_cross = response(branch A + branch B together)
              /
              [response(branch A alone) + response(branch B alone)]

recorded over all six participating sites.

This is the key attacker: the cross-branch null already contains the full
within-branch nonlinear response of each branch. I_cross therefore measures
only the additional interaction created when two nonlinear branch units are
active together.

## Modularity ruler

For each pair:

    M = mean(|log I_within,A|, |log I_within,B|)
        - |log I_cross|

Positive M means interaction is stronger **inside** the branch compartments
than **between** them.

The five-percent ruler is reused:

    M >= log(1.05)

The HUMAN primary result requires median M above that threshold, positive M in
at least 80% of all 15 pairs, and a >5% within-branch interaction in at least
four of six branches.

## NMDA attacker

The full assay is repeated with Gate 18's gamma=0.062 control rest-matched at
each site's actual settled voltage.

## Why this matters

A positive result would finally justify a stronger scaffold abstraction than
`one complicated cable`:

    many semi-independent nonlinear subunits
             -> coupled through one physical tree

That is exactly the kind of intermediate layer we need before asking whether
growth or development could construct useful scaffolds.

## Stopping line

A positive modularity result earns a branch-subunit computation assay with the
same global input budget distributed differently across compartments. It still
does not earn growth by itself.
