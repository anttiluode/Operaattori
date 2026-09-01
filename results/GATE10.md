# Gate 10 — put the matrix scaffold on a real human neuron

Gate 9 showed that ordered local operators can create effective directions that
are absent from the primitive operator set. Gate 10 deliberately does **not**
add learning yet. It asks whether the geometric half of the proposal can be
made exact on a real neuron.

## The source object

The morphology is the human L2/3 pyramidal neuron used as the high-FCI exemplar
in Aizenbud et al. (2026):

```text
morphology identifier: 1125
2013_03_06_cell11_1125_H41_06.asc
reported FCI: 0.4294
```

The file is downloaded from the authors' public `ido4848/FCI` repository at a
**pinned commit**, not from a moving branch:

```text
75ad8b4d81a7f51bf888b30650c543592340db06
```

A local file can be supplied instead, so the experiment does not require a
network connection once the morphology has been downloaded.

## Representation

The neurite morphology is flattened into a rooted point tree. Every non-root
point receives a parent-local rigid transform

```text
L_i in SE(3)
```

and world pose is recovered only by composition:

```text
W_i = W_parent(i) L_i
```

At a bifurcation one parent pose feeds multiple child transforms. The matrix
scaffold therefore splits exactly where the physical arbor splits.

The original absolute neurite coordinates are discarded from the scaffold.
Only the soma/root position, topology, local 4x4 matrices, radii and section
metadata remain.

## Gate A — exact reconstruction

Rebuild the complete morphology from the local matrices and require:

```text
max point error       < 1e-7 um
max edge-length error < 1e-7 um
rotation orthogonality/determinant error < 1e-10
```

This is intentionally a representation audit, not a learning result.

## Gate B — one matrix bends a real subtree

Choose a real internal bifurcation with a sizeable distal arbor. Change only
the rotation block of that one local matrix by 20 degrees.

The attachment point must stay fixed. Nodes outside the subtree must stay
fixed. Distal points must move. Every cable segment must retain its length.

That establishes the useful mechanical property:

> **a local matrix change propagates geometrically through the descendants
> because their coordinates live in the transported frame.**

No downstream point is individually edited.

## Why this is not Aizenbud's TCN

Aizenbud et al. use a fixed reconstructed neuron to generate I/O data and train
a fixed-size temporal convolutional network as a *complexity ruler*. Gate 10
instead changes the representation of the neuron itself: the morphology is a
branching field of transported local frames.

No claim is made that these frames are literal microtubules or a faithful
developmental model.

## Run

```bash
python -m pip install -r requirements.txt
python -m pip install -r requirements-morphology.txt
python experiments/gate10_real_neuron_scaffold.py
```

On a machine with the source ASC already downloaded:

```bash
python experiments/gate10_real_neuron_scaffold.py \
  --morphology C:/path/to/2013_03_06_cell11_1125_H41_06.asc
```

Outputs:

```text
results/gate10/gate10.json
results/gate10/cell1125_matrix_scaffold.npz
results/gate10/cell1125_original_scaffold.png
results/gate10/cell1125_one_local_bend.png
```

The NPZ is the important artifact. It contains the matrix scaffold and no
absolute neurite coordinates.

## Stopping line

A pass earns **geometry-as-scaffold** only.

Gate 11 is allowed to attach a small dynamical operator `A_i(t)` to each local
frame and ask the first functional question:

> with topology, parameter budget and local spectra controlled, does preserving
> the real spatial ordering of noncommuting local operators change a temporal
> computation relative to shuffled, commuting and dense-matrix attackers?

Only if that survives do the local frames get activity-dependent growth.
