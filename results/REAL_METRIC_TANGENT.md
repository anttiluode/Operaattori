# Real morphology metric tangent — cell 1125

## Question

Can the direct morphology compiler be differentiated with respect to **intrinsic
segment metric** on the real human cell-1125 reconstruction, and does that
analytic derivative survive the local implicit NMDA loop all the way to the
soma trace?

This is the first test in this repo where geometry is not merely an input to
the operator. It is a differentiable control variable.

## Locked object

- FCI commit: `75ad8b4d81a7f51bf888b30650c543592340db06`
- cell: human L2/3 morphology `1125`
- compiled compartments: **1653**
- branch: `PassiveDendsSimpleSomaModel[0].apic[100]`
- nonlinear sites: x = 0.25, 0.50, 0.75
- parameterization: `theta = log(local metric scale)`
- timing: 0 / 5 / 10 ms
- causal window: 40 ms at 0.025 ms
- conductance program: released HUMAN_PROBE template
- topology fixed

For one compartment at a time:

~~~text
length:
    membrane area / leak / capacitance  ~ q
    half-axial conductance              ~ 1/q

diameter:
    membrane area / leak / capacitance  ~ q
    half-axial conductance              ~ q^2

pose:
    intrinsic metric unchanged
~~~

The analytic cable derivative includes the zero-capacitance junction Schur
complement used by the direct graph compiler. The nonlinear tangent then
propagates through the same implicit Jacobian used by the causal Newton solve.

## Result

~~~text
base G reassembly relative error             1.822e-16
base C reassembly relative error             0

max dG tangent vs centered recompile          1.373e-10
max dC tangent vs centered recompile          6.464e-11

max nonlinear soma-trace tangent error        1.372e-07

Newton convergence                            yes
max Newton iterations                         2
max site consistency                          2.289e-11 mV
max tangent site consistency                  8.882e-16 mV / log-scale

pose tangent                                  exact zero
~~~

Classification:

`REAL_MORPHOLOGY_METRIC_TANGENT_VALID`

## Peak sensitivities

The baseline forward-timing soma peak was **2.355717 mV** at **27.05 ms**.

Analytic derivatives at that same peak:

| site x | local parameter | d peak / d log(q) | linearized 1% change |
|---:|---|---:|---:|
| 0.25 | length | -0.0145985 mV | -0.0001460 mV |
| 0.25 | diameter | -0.0682277 mV | -0.0006823 mV |
| 0.50 | length | -0.0382156 mV | -0.0003822 mV |
| 0.50 | diameter | -0.0381659 mV | -0.0003817 mV |
| 0.75 | length | -0.0432519 mV | -0.0004325 mV |
| 0.75 | diameter | -0.0341092 mV | -0.0003411 mV |

All six local intrinsic-metric directions are negative at this operating point:
locally adding membrane/cable load reduces the soma peak. The effect is small
for a 1% perturbation and this is **not** yet a literal growth-cone experiment.

## What this earns

The real morphology pipeline can now be written as

~~~text
morphology metric theta
        |
        v
      G(theta), C(theta)
        |
        v
 compiled passive operator
        |
        v
 local implicit NMDA
        |
        v
     soma(t; theta)

and

d soma(t; theta) / d theta
~~~

without rerunning NEURON for the metric perturbation and without
finite-differencing the nonlinear solver.

This is more specific than saying the reduced circuit is differentiable. The
**real loaded human morphology compiler itself** now has a verified intrinsic
metric tangent.

## Fence

This validates differentiation of the compiled cable abstraction.

It does **not** show that biological dendrites optimize soma peak, that growth
uses gradient descent, or that intracellular organelles can be ignored in
general. Intracellular variables belong in this abstraction if they change the
effective membrane, axial, channel, synaptic, or topology state.

The next useful move is a sensitivity map over many compartments using this
analytic tangent, followed by a literal distal-growth/load perturbation rather
than another generic geometry claim.
