# Gate 24 — real extracellular field x movable scaffold x nonlinear branch

Gate 23 was an abstract world-sampling toy. Gate 24 is the first scaffold gate
that applies **extracellular potential through NEURON's cable mechanism**.

The broad biological mechanism is not claimed as new.

- Aspart, Remme & Obermayer (2018) derive morphology-dependent polarization of
  straight and bent cables in extracellular fields.
- Fan et al. (2023/2024) explicitly model electric-field regulation of NMDA
  dendritic integration and NMDA-spike generation.

Operaattori asks a narrower question:

> Holding the intrinsic cell-1125 cable exactly fixed, can an isometric local
> matrix bend change the nonlinear interaction of one of the branch
> compartments already earned in Gates 16–20?

## Released model

Use the pinned FCI cell-1125 NEURON model from Gates 16–21.

The dendrites in this released model are passive. Therefore Gate 24 does **not**
claim to reproduce the strong active 10–20 Hz apical resonance that Aspart et
al. attribute largely to h-type conductance. We nevertheless include 10, 15
and 20 Hz as a locked low-frequency panel.

## Geometry

Recover the same six compact three-site branch sections used by Gate 20.

For each candidate section, compute its complete descendant-section cable
length before any field response is seen.

Two deterministic arms are then chosen:

- **PROXIMAL / LARGE-SUBTREE** — largest descendant cable extent;
- **DISTAL / SMALL-SUBTREE** — smallest descendant cable extent.

For either arm, the bend is a rigid 35-degree rotation of that section and all
of its descendants about the section's proximal attachment point.

NEURON's section lengths, diameters, topology and mechanisms are never changed.
Only the XYZ coordinates used to compute extracellular potential are rotated.

This is the Gate-22 isometry idea applied directly to the released dynamic
model.

## Field

Apply one spatially uniform field of amplitude

    E0 = 1 V/m

along the morphology's deterministic apical principal axis.

For a material point at world coordinate r,

    Vext(r,t) = - E(t) d dot (r - r_soma)

with the usual unit conversion to mV.

The AC waveform is phase-locked so E(t) is at its positive peak at the synaptic
event. Frequencies are locked to:

    10, 15, 20 Hz

Unlike Gate 23, extracellular potential is **not normalized by total sampled
drive**. A uniform electric field is the physical control; renormalizing its
spatial profile would change the experiment.

## Nonlinear readout

Use Gate 20's locked compact dose:

    3 sites x multiplicity 8 = 24 simultaneous virtual synapses

For each branch, geometry, frequency and NMDA condition:

1. run the field with no synaptic event;
2. run all three sites together;
3. run each of the three sites alone;
4. subtract the field-only voltage trace from every synaptic run;
5. compare the simultaneous branch response with the sum of the three
   single-site perturbations.

Define

    I = AUC(simultaneous local synaptic perturbation)
        ------------------------------------------------
        AUC(sum of exact single-site perturbations)

and the bend modulation

    B = |log(I_bent / I_original)|.

The same assay is run for HUMAN gamma=0.078 and Gate 20's sitewise
rest-matched gamma=0.062 attacker.

## Controls

### Zero field

With drive set to zero, ORIGINAL and BENT must be identical. The cable has not
changed.

### Material locked

Keep the ORIGINAL extracellular coefficients attached to the same material
segments while calling the geometry "bent." This must also be identical.

### Distal arm

The small-subtree bend is the anatomical attacker against the large-subtree
arm.

### Field-only trace

All nonlinear measurements are made after subtracting the exact field-only
trace. A field-induced baseline oscillation cannot masquerade as synaptic
interaction.

## Preregistered primary prediction

Before running:

1. zero-field and material-locked relative differences must be < 1e-6;
2. the large-subtree HUMAN bend must change the nonlinear interaction by at
   least 5% (B >= log(1.05)) at two of the three frequencies;
3. its median B must be at least 2x the small-subtree arm's median B;
4. no somatic spike may contaminate the branch-interaction assay.

If all four hold:

    REAL_FIELD_BEND_MODULATES_NONLINEAR_COMPARTMENT

If the extracellular field polarizes the branch but the nonlinear interaction
does not clear the ruler:

    FIELD_COUPLING_PRESENT_BUT_NOT_COMPARTMENT_SELECTIVE

Otherwise:

    NO_ROBUST_REAL_FIELD_BEND_EFFECT

The rest-matched gamma arm is reported as a mechanism attacker but is not
required to be smaller for the primary scaffold claim. Prior literature already
establishes that electric fields can regulate NMDA integration; this gate is
about **isometric scaffold position**.

## Literature fence

Aspart et al. 2018:
https://doi.org/10.1371/journal.pcbi.1006124

Fan et al. 2023/2024:
https://doi.org/10.1007/s11571-022-09922-y

NEURON extracellular mechanism:
https://nrn.readthedocs.io/

## Stopping line

No Gate 25 is opened in this session.

If Gate 24 passes, the result is demo-worthy but still not a growth result.
If it fails, do not increase field amplitude, scan bend angles or tune synaptic
dose to force a threshold crossing.


## Receipt — real field couples, compartment selectivity does not

The corrected locked CI assay completed successfully.

~~~text
distal small-subtree arm
  section                              apic[58]
  descendant cable                     316.9 um
  max moved segment                    159.9 um

proximal large-subtree arm
  section                              apic[77]
  descendant cable                     599.3 um
  max moved segment                    278.6 um

controls
  zero-field max relative difference   0
  material-locked max difference       0

HUMAN nonlinear branch interaction
  proximal median bend factor          1.0003x
  distal median bend factor            1.0002x
  proximal frequencies >5%             0 / 3
  proximal/distal effect ratio         1.206x
  spike guard fraction                 0
~~~

Classification:

~~~text
FIELD_COUPLING_PRESENT_BUT_NOT_COMPARTMENT_SELECTIVE
~~~

The causal wiring passed its strongest bookkeeping checks exactly: when the
field is removed, changing only the embedding is invisible; when the original
field coefficients are locked to material segments, the nominal bend is also
invisible.

With a real uniform 1 V/m extracellular field, however, the locked 35-degree
re-embedding changes the already-measured compact-branch nonlinear interaction
by only about three parts in ten thousand in the large-subtree arm. No tested
frequency clears the preregistered 5% ruler.

So Gate 24 completes rather than rescues the scaffold decomposition:

> real extracellular field coupling exists, but this weak-field isometric bend
> does not selectively recruit the nonlinear branch compartment under the
> locked protocol.

Do not increase field amplitude, scan bend angle, or tune synaptic dose after
seeing this result. Gate 25 remains closed.
