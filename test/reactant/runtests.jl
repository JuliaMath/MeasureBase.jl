# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Reactant smoke tests, not part of the default test suite. Run with
# `julia --project=test/reactant test/reactant/runtests.jl` after
# instantiating that project, or include this file in an environment that
# provides Reactant.

using Test
using Reactant
using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, StdLogistic, Lebesgue, Dirac
using MeasureBase: logdensities, logdensity_rel, weightedmeasure, superpose, restrict, mintegrate_exp
using MeasureBase: mcombine
using ArraysOfArrays: VectorOfSimilarVectors, sliced, flatview
using Distributions: Normal, Exponential, Uniform, Beta

Reactant.set_default_backend("cpu")

# Compiles `f` for traced copies of `args` and compares with the plain result:
function test_traced(f, args...; kwargs...)
    expected = f(args...)
    traced_args = map(Reactant.to_rarray, args)
    result = @jit f(traced_args...)
    @test _plain(result) ≈ _plain(expected) nans = true
    return result
end

_plain(x::AbstractArray) = Array(x)
_plain(x::Number) = Float64(x)

@testset "Reactant" begin
    x = randn(10)
    X = randn(3, 20)

    @testset "powers and batches" begin
        test_traced(x -> logdensityof(StdNormal()^10, x), x)
        test_traced(X -> logdensities(StdNormal(), X), X)
        test_traced(X -> logdensities(StdNormal()^3, X), X)
        test_traced(X -> logdensities(StdNormal()^3, sliced(X, 1)), X)
        test_traced(X -> logdensities((StdNormal()^3)^4, reshape(X[:, 1:16], 3, 4, 4)), X)
        test_traced(x -> logdensityof(StdUniform()^10, x), rand(10))
        test_traced(x -> logdensityof(StdExponential()^10, x), rand(10))
        test_traced(x -> logdensityof(weightedmeasure(0.3, StdNormal())^10, x), x)
        test_traced(x -> logdensityof(Lebesgue()^10, x), x)
    end

    @testset "support masks" begin
        xu = 2 .* rand(10) .- 0.5
        test_traced(x -> logdensities(StdUniform(), x), xu)
        test_traced(x -> logdensities(StdExponential(), x), xu)
        test_traced(x -> logdensities(restrict(x -> x > 0, StdNormal()), x), xu)
    end

    @testset "relative densities" begin
        xu = 2 .* rand(10) .- 0.5
        test_traced(x -> logdensity_rel.(Ref(StdUniform()), Ref(StdExponential()), x), xu)
        test_traced(x -> logdensity_rel.(Ref(StdNormal()), Ref(StdLogistic()), x), x)
        test_traced(x -> logdensity_rel.(Ref(StdNormal()^10), Ref(StdLogistic()^10), Ref(x)), x)
    end

    @testset "superposition, density measures and spike mixtures" begin
        mix = superpose(weightedmeasure(log(0.3), StdNormal()), weightedmeasure(log(0.7), StdLogistic()))
        test_traced(x -> logdensities(mix, x), x)
        test_traced(x -> logdensityof(mix^10, x), x)
        dm = mintegrate_exp(x -> -abs(x), StdNormal())
        test_traced(x -> logdensities(dm, x), x)
        sm = SpikeMixture(StdNormal(), 0.2)
        test_traced(x -> logdensities(sm, x), vcat(x, 0.0))
    end

    # Products over arrays of marginals are not covered: Reactant can't
    # broadcast over arrays of measures together with traced arrays.
    @testset "structural batched kernels" begin
        w = weightedmeasure(log(0.3), StdNormal()^3)
        test_traced(X -> logdensities(w, X), X)
        mc = mcombine(vcat, StdNormal()^2, StdUniform()^1)
        test_traced(X -> logdensities(mc, X), vcat(X[1:2, :], rand(1, 20)))
        test_traced(X -> logdensities((StdNormal()^2)^3, X), reshape(X[1:2, 1:6], 2, 3, 2))
    end

    @testset "transport of powers and products" begin
        test_traced(x -> transport_to(StdNormal()^10, StdUniform()^10)(x), rand(10))
        test_traced(x -> transport_to(StdExponential()^10, StdNormal()^10)(x), x)
        test_traced(x -> transport_to(StdLogistic()^10, StdExponential()^10)(x), rand(10))
        test_traced(x -> transport_to(StdNormal()^6, (StdUniform()^2)^3)(x), rand(2, 3))
        mc = mcombine(vcat, StdNormal()^2, StdUniform()^1)
        test_traced(x -> transport_to(StdNormal()^3, mc)(x), vcat(randn(2), rand(1)))
        test_traced(z -> transport_to(mc, StdNormal()^3)(z), randn(3))
    end

    @testset "batched transport" begin
        test_traced(X -> transport_to(StdNormal(), StdUniform()).(X), rand(10))
        test_traced(X -> flatview(transport_to(StdExponential()^3, StdNormal()^3).(sliced(X, Val(1)))), X)
        mc = mcombine(vcat, StdNormal()^2, StdUniform()^1)
        test_traced(X -> flatview(transport_to(StdLogistic()^3, mc).(sliced(X, Val(1)))), vcat(X[1:2, :], rand(1, 20)))
    end

    @testset "transport" begin
        test_traced(x -> transport_to(StdUniform(), StdNormal()).(x), x)
        test_traced(x -> transport_to(StdNormal(), StdUniform()).(x), rand(10))
        test_traced(x -> transport_to(Normal(2, 3), StdNormal()).(x), x)
        test_traced(x -> transport_to(StdNormal(), Exponential(2.0)).(x), rand(10))
    end
end
