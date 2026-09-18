# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Reactant tests. Reactant isn't a static test dependency (it only
# supports 64-bit Linux and macOS), runtests.jl adds it on the fly where
# supported. The backend defaults to the CPU, set the environment variable
# `MEASUREBASE_REACTANT_BACKEND` (e.g. to "gpu") to change it; the file can
# also be run standalone in an environment that provides Reactant.

using Test
using Reactant
using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, StdLogistic, Lebesgue, Dirac, asmeasure
using MeasureBase: logdensities, logdensity_rel, weightedmeasure, superpose, restrict, mintegrate_exp
using MeasureBase: mcombine
using ArraysOfArrays: VectorOfSimilarVectors, sliced, flatview
using Distributions: Normal, Uniform, Exponential, Logistic, Cauchy, Laplace, LogNormal, Weibull, Gamma, Beta
using Distributions: Poisson, MvNormal, Dirichlet
using MeasureBase.InverseFunctions: inverse

Reactant.set_default_backend(get(ENV, "MEASUREBASE_REACTANT_BACKEND", "cpu"))

# Compiles `f` for traced copies of `args` and compares with the plain
# result. Array results are copied inside the compiled function, so that
# views and reshapes of device arrays come back as plain device arrays:
function test_traced(f, args...; kwargs...)
    expected = f(args...)
    traced_args = map(Reactant.to_rarray, args)
    g = (xs...) -> _contiguous(f(xs...))
    result = @jit g(traced_args...)
    @test _plain(result) ≈ _plain(expected) nans = true
    return result
end

_contiguous(x::AbstractArray) = copy(x)
_contiguous(x) = x
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

    # Distribution parameters stay constants, Distributions' parameter
    # structs can't hold traced arrays:
    @testset "wrapped distributions" begin
        for d in (Normal(0.3, 1.7), Uniform(-1.0, 2.5), Exponential(0.7), Logistic(0.2, 1.3), Cauchy(0.1, 0.8), Laplace(-0.4, 1.1), LogNormal(0.2, 0.6), Weibull(1.4, 0.9), Gamma(2.3, 1.2), Beta(2.5, 3.5))
            m = asmeasure(d)
            xd = rand(d, 10)
            test_traced(x -> logdensities(m, x), xd)
            f = transport_to(StdNormal(), m)
            if d isa Union{Gamma,Beta}
                # SpecialFunctions' incomplete gamma and beta functions have no Reactant methods:
                @test_broken @jit((x -> copy(f.(x)))(Reactant.to_rarray(xd))) isa AbstractArray
            else
                test_traced(x -> f.(x), xd)
                test_traced(z -> inverse(f).(z), randn(10))
            end
        end
        test_traced(x -> logdensities(asmeasure(Poisson(2.7)), x), Float64.(rand(Poisson(2.7), 10)))
        mvn = MvNormal([0.3, -2.9], [1.7 0.5; 0.5 2.3])
        mm = asmeasure(mvn)
        Xm = rand(mvn, 10)
        test_traced(X -> logdensities(mm, X), Xm)
        test_traced(X -> flatview(transport_to(StdNormal()^2, mm).(sliced(X, Val(1)))), Xm)
        test_traced(Z -> flatview(transport_to(mm, StdNormal()^2).(sliced(Z, Val(1)))), randn(2, 10))
        dir = Dirichlet([2.0, 3.0, 4.0, 1.5])
        md = asmeasure(dir)
        test_traced(X -> logdensities(md, X), rand(dir, 10))
    end
end
