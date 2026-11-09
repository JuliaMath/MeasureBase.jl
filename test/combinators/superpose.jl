using Test

using MeasureBase
using MeasureBase: superpose, weightedmeasure, StdNormal, StdUniform, SpikeMixture, GenContext
using Distributions: Poisson
using MeasureBase.HeterogeneousComputing: AbstractComputeUnit
using JLArrays: JLArray, JLBackend
using Random: Xoshiro

@testset "superpose.jl" begin
    μ = Dirac(0)
    ν = Dirac(1)
    μs = μ + ν
    @test μs isa SuperpositionMeasure{<:Tuple{Dirac,Dirac}}
    @test μs == SuperpositionMeasure((μ, ν)) == superpose(μ, ν)
    @test exp(logdensity_def(μs, 0)) == 1.0
    @test basemeasure(μs) == μs
    @test densityof(μs, 0) == 1.0
    @test densityof(μs, 1) == 1.0
    @test densityof(μs, 0.5) == 0

    μs = SuperpositionMeasure([μ, ν])
    @test μs isa SuperpositionMeasure{<:AbstractVector{<:AbstractMeasure}}
    @test exp(logdensity_def(μs, 0)) == 1.0
    @test collect(basemeasure(μs).components) == [μ, ν]
    @test densityof(μs, 0) == 1.0
    # Base measures of components count wherever they have mass, not only
    # where the component itself does:
    @test logdensityof(superpose(StdNormal(), StdUniform()), -1.0) ≈ logdensityof(StdNormal(), -1.0)
    @test logdensityof(superpose(StdNormal(), StdUniform()), 0.5) ≈ log(exp(logdensityof(StdNormal(), 0.5)) + 1)

    # Discrete components: bases carry their atoms, densities are relative
    # to the superposition of Lebesgue and counting measures:
    N = StdNormal()
    ld_N(x) = logdensityof(N, x)
    @test logdensityof(superpose(N, StdUniform()), 1.5) ≈ ld_N(1.5)
    sd = N + 0.3 * Dirac(0.0)
    @test logdensityof(sd, 0.0) ≈ log(0.3)
    @test logdensityof(sd, 0.5) ≈ ld_N(0.5)
    @test logdensities(sd, [0.0, 0.5]) ≈ [log(0.3), ld_N(0.5)]
    sc = N + Counting()
    @test logdensityof(sc, 1.0) ≈ 0
    @test logdensityof(sc, 0.5) ≈ ld_N(0.5)
    @test logdensity_rel(Dirac(0.0) + Counting(), Counting(), 0.0) ≈ log(2)
    @test logdensity_rel(Dirac(0.0) + Counting(), Counting(), 1.0) ≈ 0
    @test logdensityof(Dirac(0.0) + Counting(), 0.5) == -Inf
    # Spike and slab:
    ss = 0.4 * N + 0.6 * Dirac(0.0)
    sm = SpikeMixture(N, 0.4)
    for x in (0.0, 0.7, -1.2)
        @test logdensityof(ss, x) ≈ logdensityof(sm, x)
    end

    # Spike mixtures with spike weights independent of the parent weight:
    sm2 = SpikeMixture(N, 0.3, 0.5)
    @test massof(sm2) ≈ 0.8
    @test logdensityof(sm2, 0.0) ≈ log(0.5)
    @test logdensityof(sm2, 1.0) ≈ log(0.3) + ld_N(1.0)
    @test logdensities(sm2, [0.0, 1.0]) ≈ [log(0.5), log(0.3) + ld_N(1.0)]
    draws = [rand(sm2) for _ in 1:10^5]
    @test isapprox(count(iszero, draws) / length(draws), 0.5 / 0.8, atol = 0.01)
    @test isapprox(count(iszero, rand(sm2^(10^5))) / 10^5, 0.5 / 0.8, atol = 0.01)

    @test MeasureBase.insupportof(MeasureBase.TrivialMeasure(), 1.0) == false
    @test logdensity_rel(MeasureBase.TrivialMeasure(), MeasureBase.LebesgueBase(), 1.0) == -Inf

    # Dirac equality is not decidable from types, so no weighted collapse:
    μ2 = μ + μ
    @test μ2 isa SuperpositionMeasure
    @test μ2 == superpose(μ, μ)
    @test exp(logdensity_def(μ2, 0)) == 1.0

    # Relative densities of superpositions of equal type don't assume equal
    # components, the oracle is the mass of the atoms:
    let a = superpose((Dirac(0.0), Dirac(0.0))), b = superpose((Dirac(0.0), Dirac(1.0)))
        @test logdensity_rel(a, b, 0.0) ≈ log(2)
        @test logdensity_rel(b, a, 0.0) ≈ -log(2)
        @test logdensity_rel(b, b, 0.0) == 0
        @test logdensity_rel(b, a, 1.0) == Inf
        @test logdensity_rel(a, b, 1.0) == -Inf
    end
    let a = SuperpositionMeasure((2 * StdNormal(), 3 * StdNormal()))
        b = SuperpositionMeasure((StdNormal(), StdNormal()))
        @test logdensity_rel(a, b, 0.3) ≈ log(5 / 2)
    end

    # For singleton measure types equal measures combine into weighted measures:
    s2 = StdNormal() + StdNormal()
    @test s2 isa WeightedMeasure
    @test exp(s2.logweight) ≈ 2
    @test basemeasure(s2) == StdNormal()
end

@testset "random variates of superpositions" begin
    # Draws of the components are promoted to a common type, independent of
    # which components are selected:
    s = superpose(StdNormal(), Dirac(0.0))
    for seed in 1:6
        @test MeasureBase.batched_rand_impl(GenContext{Float32}(Xoshiro(seed)), s, (4,)) isa Vector{Float64}
    end
    @test @inferred(MeasureBase.batched_rand_impl(GenContext{Float32}(Xoshiro(1)), s, (4,))) isa Vector{Float64}
    @test @inferred(rand(Xoshiro(1), s)) isa Float64
    @test rand(superpose(StdNormal(), Dirac(0))^4) isa Vector{Float64}
    @test eltype(rand(superpose(asmeasure(Poisson(2.0)), Dirac(0.5))^4)) == Float64
    dctx = GenContext{Float32}(AbstractComputeUnit(JLBackend()))
    @test MeasureBase.batched_rand_impl(dctx, s, (5,)) isa JLArray{Float64}

    # Components are selected by mass, single variates like batches:
    m = superpose(0.3 * StdNormal(), 0.7 * Dirac(2.0))
    X = rand(Xoshiro(1), m^100_000)
    @test count(==(2.0), X) / length(X) ≈ 0.7 atol = 0.01
    x = [rand(Xoshiro(i), m) for i in 1:10_000]
    @test count(==(2.0), x) / length(x) ≈ 0.7 atol = 0.02
    v = rand(Xoshiro(1), superpose(StdNormal()^2, Dirac([0.0, 0.0])))
    @test v isa Vector{Float64} && length(v) == 2
end
