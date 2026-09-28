# MeasureBase.jl

This package contains core definitions for (probability and non-probability) measures.

MeasureBase provides primitive and standard measures, measure combinators
(powers, products, pushforwards, binds, superpositions and more),
log-densities, transports between measures and random variate generation.
Everything works for single variates and for batches of variates, on the
CPU, on GPUs and under Reactant.

See [Concepts](@ref) for an overview of the design, the conventions and
the extension points, and the [API](@ref) for the reference documentation.
