# Gate 13 — put branch junctions into the transport algebra\n\nGate 12 killed the fine serial-order interpretation: the real root-to-tip\npaths are already ordinary thick-to-thin taper.\n\nThe next neuron-specific operation is therefore not another path shuffle. It\nis the bifurcation.\n\n## Branch algebra\n\nA passive cable segment is represented by its exact two-port matrix. At a\nbranch point, daughter subtrees contribute input admittances in parallel.\nAlong one selected soma-to-tip route, every off-path subtree is therefore an\nexact shunt operator inserted between cable operators:\n\n    M1 S(Yside,1) M2 S(Yside,2) ... Mn\n\nwhere\n\n    S(Y) = [[1, 0],\n            [Y, 1]]\n\nGate 13 computes the same transfer two ways:\n\n1. exact whole-tree Schur elimination;\n2. explicit cable-plus-side-shunt matrix composition along each tip path.\n\nThose must agree before any branch effect is interpreted.\n\n## Causal ablation\n\nFor every dendritic tip, compare:\n\n- FULL: the real tree, including all off-path daughter loads;\n- PATH ONLY: the exact same soma-to-tip cable pieces with every side branch removed.\n\nThe path taper is therefore held fixed. Only branch loading is removed.\n\nRead gain, phase, complex transfer difference, and the effective rank of the\npopulation of tip-to-soma transfer signatures.\n\n## What this gate can earn\n\nA positive result means branch junctions contribute a functional operator\nbeyond serial tapering paths. It does not establish nonlinear dendritic\ncomputation, learning, or biological optimality.\n\nIf branch loading is substantial, the natural next gate is to add the\nnonlinear synaptic ingredient emphasized by Aizenbud rather than inventing an\narbitrary neural network inside the scaffold.\n

## First receipt — branch loading is real, passive portfolio remains low-rank

The exact full-tree job passed. Whole-tree Schur elimination and the
independent cable-plus-side-shunt matrix product agree to numerical precision:

    dendritic tips                         110
    median path length                     281.8 um
    median branch junctions/path           5.0
    tree vs shunt-product max rel error    4.558e-11

    full vs isolated median difference     0.2676
    median |gain change|                   1.712 dB
    median signed gain change             -1.712 dB
    median phase change                    0.1660 rad
    tip-frequency points >10% changed      0.736
    tips with median effect >10%            0.864
    branch-count/effect correlation        0.958

    signature effective rank FULL          1.030
    signature effective rank PATH          1.024

Classification:

    BRANCH_LOADING_CHANGES_TRANSFER

The causal ablation keeps each tested soma-to-tip tapering path identical and
only removes off-path daughter subtrees. So unlike Gate 11, this result cannot
be reduced to reversing the gross taper. Branch junctions themselves are a
functional part of the passive operator scaffold.

The negative result is equally important: the population of passive tip
transfer signatures remains almost rank one. Branch loading changes amplitudes
and phases strongly but does not by itself create a high-dimensional repertoire
of independent computations.
