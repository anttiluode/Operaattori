# Operaattori prior-art kill ledger — September 2026

This file asks a deliberately hostile question:

> Which Operaattori claims remain scientifically distinctive after comparing them with the strongest nearby literature?

It is a claim-boundary document, not a novelty declaration. `NOT KILLED IN THIS PASS` means only that this search did not locate the same combination; it does **not** mean priority has been established.

## Bottom line

Most broad interpretations of Operaattori are **not novel**. The defensible object, if any, is narrow:

> an explicit zero-fit morphology-to-cable compiler whose first-order intrinsic-metric derivative is propagated analytically through junction elimination and a local implicit voltage-dependent NMDA closure, with a direct transport-versus-NMDA-feedback decomposition measured across the released 24-cell FCI morphology panel.

Even that should be presented as a methods / quantitative decomposition contribution unless a deeper literature search establishes stronger priority.

## Kill table

| Operaattori statement | Strong nearby prior art | Verdict |
|---|---|---|
| Dendritic morphology compiles a passive transfer operator | Rall cable theory; Koch & Poggio, *J Neurosci Methods* 1985, algorithms for arbitrary dendritic trees; decades of transfer-impedance work | **KILLED as conceptual novelty** |
| A morphology can be reduced while preserving dendro-somatic transfer and nonlinear dendritic behavior | Amsalem et al., **Neuron_Reduce**, *Nature Communications* 2020, doi:10.1038/s41467-019-13932-6 | **KILLED as broad reduction claim** |
| A small set of Green/transfer kernels can stand in for a full morphological cable when only a few sites matter | Wybo, Stiefel & Torben-Nielsen, *Biological Cybernetics* 2013, doi:10.1007/s00422-013-0568-0; Wybo et al., *Neural Computation* 2015, doi:10.1162/NECO_a_00788 | **KILLED as Green-kernel concept** |
| Detailed biophysical neurons and their morphology can be differentiated / optimized | Deistler et al., **Jaxley**, *Nature Methods* 2025, doi:10.1038/s41592-025-02895-w; Jaxley explicitly optimizes compartment length and radius | **KILLED** |
| Length/diameter sensitivity can vary in magnitude and sign with operating point / active state | Weaver & Wearne, *PLoS Computational Biology* 2008, doi:10.1371/journal.pcbi.0040011; earlier and parallel neuronal sensitivity-analysis literature | **KILLED as general principle** |
| Passive dendritic transfer has geometry-dependent sensitivity structure | Korogod & Kaspirzhny, *Biological Cybernetics* 2008, doi:10.1007/s00422-007-0204-y | **KILLED** |
| Morphology changes NMDA integration because local/somatic input resistance changes | Poleg-Polsky, *PLoS ONE* 2015, doi:10.1371/journal.pone.0140254 | **KILLED** |
| Thin/distal versus thick/proximal dendrites can reverse the relative effectiveness of passive/AMPA and voltage-dependent NMDA drive | Lajeunesse, Kröger & Timofeev, *Journal of Neurophysiology* 2013, doi:10.1152/jn.01090.2011 | **KILLED as qualitative mechanism** |
| Morphology and NMDA nonlinearity jointly shape computation in the released human/rat cortical morphology family | Aizenbud et al., *PNAS* 2026, doi:10.1073/pnas.2533168123, the FCI work from which the 24-model panel derives | **KILLED as biological headline** |
| Rigid 3-D pose changes nothing if intrinsic cable metric/topology is fixed | Classical cable theory depends on intrinsic electrical geometry rather than arbitrary Euclidean embedding unless an embedding-dependent mechanism is explicitly added | **KILLED as scientific novelty; useful implementation fence only** |
| A causal local implicit solve is more numerically portable than one global waveform Picard iteration for this reduced NMDA circuit | Numerical analysis strongly anticipates state-local implicit methods outperforming fragile global fixed-point iterations near strong nonlinear feedback | **Useful engineering result, not presently a novelty claim** |
| The exact derivative can be written as a passive/transport contribution plus a voltage-dependent NMDA feedback contribution | Chain rule / implicit-function sensitivity makes such a split mathematically natural; active-versus-morphology interactions and NMDA morphology effects are old | **Broad principle killed** |
| Across the 24 released morphologies, the *measured local geometry derivative* is decomposed into frozen-current transport and NMDA-feedback terms, and feedback reverses at least one tested geometry direction in 23/24 cells | No exact matching paper was located in this pass. Prior work establishes the ingredients and qualitative mechanism, but not this exact derivative decomposition and cross-cell panel | **NOT KILLED IN THIS PASS; narrow quantitative seam** |
| The direct graph compiler plus hand-derived tangent through zero-capacitance junction Schur elimination and local implicit NMDA gives ~1e-7 tangent agreement on a 1653-compartment released human reconstruction without fitting | Jaxley provides general autodiff; Neuron_Reduce / Green-function work provide reduction; classical sensitivity work provides derivatives. No exact same transparent compiler+tangent validation was located here | **NOT KILLED IN THIS PASS; methods seam, priority unproven** |

## Why the 24-cell decomposition is narrower than the old story

The broad sentence

> morphology supplies transport and NMDA state changes what morphology does

is old territory.

The surviving measurable statement is more specific. Operaattori defines, at the same operating point and for the same local metric perturbation,

```text
full geometry derivative
    =
frozen-current transport derivative
    +
voltage-dependent feedback derivative
```

and measures the vectors separately rather than inferring their roles from two different cell classes or receptor conditions.

The current locked panel reports:

```text
24 released morphologies
3 drives
6 local metric directions per cell per drive

median feedback/transport norm:
    0.827 at 0.5x
    1.468 at 1.0x
    1.735 at 2.0x

cells with feedback > transport at >=1 drive: 24/24
cells with >=1 full-vs-transport sign reversal: 23/24
```

At 1x, 38 of the 39 sign reversals are diameter directions. This is a quantitative decomposition result, not evidence that NMDA/morphology interaction itself was unknown.

## Why the analytic tangent is narrower than Jaxley

Jaxley removes any claim that differentiable morphology or gradient-based shape optimization is new.

Operaattori's remaining methods distinction is implementation transparency rather than differentiability itself:

```text
intrinsic metric
   -> exact local dG,dC
   -> junction Schur/elimination derivative
   -> compiled passive state-space operator derivative
   -> local implicit NMDA Jacobian reuse
   -> soma-trace derivative
```

On cell 1125 (1653 compartments), the archived receipt reports worst nonlinear soma-trace tangent error about `1.37e-7` against centered metric recompilation. This is useful validation of a hand-derived reduced compiler. It should not be advertised as the first differentiable neuron simulator or first morphology optimizer.

## The strongest prior-art attacks

### Jaxley (2025)

Deistler et al. make detailed biophysical models differentiable with JAX/autodiff and explicitly train a single-neuron morphology including **length and radius**. Therefore all claims of “first differentiable dendritic morphology,” “geometry as a trainable variable,” or “gradient-based morphology optimization” are out.

### Green-function formalism (2013/2015)

Wybo and colleagues already use analytic cable Green functions to compress the influence of a detailed morphology into interaction kernels when only a limited set of locations matters. The 2015 sparse reformulation is especially close to the idea that only the selected input/output sites need participate in the reduced dynamical object.

Therefore `local Green matrix × synapse law` is not, by itself, a new conceptual architecture.

### Neuron_Reduce (2020)

Amsalem et al. analytically reduce detailed nonlinear neurons while preserving transfer impedance, remap synapses and active channels, and retain nonlinear dendritic computations. Therefore “zero/low-fit analytical reduction preserving transfer” already has strong prior art.

### Morphology sensitivity (2007–2008)

Weaver/Wearne and Korogod/Kaspirzhny already perform mathematical sensitivity analysis on dendritic morphology / transfer and show that sensitivity depends on morphology, passive state, and active conductances. Sign-changing / operating-point-dependent morphology sensitivity is not an Operaattori discovery.

### NMDA × morphology (2013/2015 and later)

Lajeunesse et al. show that AMPA and NMDA somatic efficacy follow different morphology-dependent regimes in reconstructed thalamocortical trees. Poleg-Polsky explicitly varies dendritic diameter and length and shows how input resistance changes NMDA-spike thresholds and somatic effects. These papers are close enough that the qualitative transport-versus-feedback intuition cannot carry novelty.

### FCI parent paper (2026)

Aizenbud et al. already establish that morphology and the density/nonlinearity of NMDA signaling jointly explain functional complexity differences in the same 24 human/rat cortical model collection. Operaattori may offer a mechanistic operator decomposition of those models; it cannot claim the parent biological conclusion as its own.

## Claim language that survives this audit

Reasonable:

> We built a transparent zero-fit cable compiler and an analytic intrinsic-metric tangent through a local implicit NMDA closure, then used it to decompose local morphology sensitivity into transport and voltage-dependent feedback on a released 24-morphology panel.

> In the tested panel, nonlinear feedback becomes comparable to or larger than passive transport and reverses the sign of at least one tested local geometry direction in 23/24 cells.

Not reasonable:

> We discovered that dendritic geometry matters.

> We discovered that NMDA nonlinearity interacts with morphology.

> We discovered differentiable dendritic morphology.

> We discovered a Green-function reduction of neurons.

> We discovered that a neuron's 3-D pose is irrelevant to classical cable transfer.

## What would turn the remaining seam into a real contribution?

Do **not** add another conceptual gate. Compare the decomposition directly against the closest established tools.

1. **Jaxley equivalence attack.** Recreate a locked subset of the 24-cell metric derivatives in Jaxley/autodiff. The question is not whether both gradients exist, but whether Operaattori's explicit transport/feedback split is numerically reproducible from a general differentiable simulator.
2. **Green-function / reduction attack.** Demonstrate exactly what information Operaattori retains that a standard Green-function or transfer-impedance reduction does not when NMDA feedback is strong.
3. **Prediction attack.** Use the decomposition *before* a held-out finite metric perturbation to predict where the passive/transport sign will be reversed by nonlinear feedback. A decomposition that only explains an already-computed gradient is less interesting than one that predicts a finite intervention.
4. **Cross-model attack.** Repeat the split in at least one non-FCI detailed model family. This separates a reusable method from an audit specialized to one released model set.
5. **Cost/accuracy attack.** Benchmark runtime, memory, gradient accuracy, and numerical robustness against Jaxley and/or NEURON finite differences at matched observables.

If Operaattori cannot win or uniquely diagnose anything on those comparisons, classify it as a useful transparent educational/research implementation rather than a new method.

## References

- Deistler M. et al. **Jaxley: differentiable simulation enables large-scale training of detailed biophysical models of neural dynamics.** *Nature Methods* 22, 2649–2657 (2025). doi:10.1038/s41592-025-02895-w
- Amsalem O. et al. **An efficient analytical reduction of detailed nonlinear neuron models.** *Nature Communications* 11, 288 (2020). doi:10.1038/s41467-019-13932-6
- Wybo W.A.M., Stiefel K.M., Torben-Nielsen B. **The Green's function formalism as a bridge between single- and multi-compartmental modeling.** *Biological Cybernetics* 107, 685–694 (2013). doi:10.1007/s00422-013-0568-0
- Wybo W.A.M. et al. **A Sparse Reformulation of the Green's Function Formalism Allows Efficient Simulations of Morphological Neuron Models.** *Neural Computation* 27, 2587–2622 (2015). doi:10.1162/NECO_a_00788
- Weaver C.M., Wearne S.L. **Neuronal Firing Sensitivity to Morphologic and Active Membrane Parameters.** *PLoS Computational Biology* 4, e11 (2008). doi:10.1371/journal.pcbi.0040011
- Korogod S.M., Kaspirzhny A.V. **Parameter sensitivity of distributed transfer properties of neuronal dendrites: a passive cable approximation.** *Biological Cybernetics* 98, 87–100 (2008). doi:10.1007/s00422-007-0204-y
- Lajeunesse F., Kröger H., Timofeev I. **Regulation of AMPA and NMDA receptor-mediated EPSPs in dendritic trees of thalamocortical cells.** *Journal of Neurophysiology* 109, 13–30 (2013). doi:10.1152/jn.01090.2011
- Poleg-Polsky A. **Effects of Neural Morphology and Input Distribution on Synaptic Processing by Global and Focal NMDA-Spikes.** *PLoS ONE* 10, e0140254 (2015). doi:10.1371/journal.pone.0140254
- Aizenbud I. et al. **Dendritic morphology and synaptic nonlinearities enhance functional complexity in human cortical neurons.** *PNAS* 123, e2533168123 (2026). doi:10.1073/pnas.2533168123
