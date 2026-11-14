using Test

using MeasureBase
using MeasureBase: superpose, weightedmeasure, StdNormal, StdUniform, StdExponential, GenContext
using MeasureBase: MeasureJoin, LebesgueBase, Lebesgue, static
using Distributions: Distributions, Poisson, MvNormal, Dirichlet, Normal, pdf, logpdf
using StaticArrays: SVector
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
    @test basemeasure(μs) == MeasureJoin((μ, ν))
    @test densityof(μs, 0) == 1.0
    @test densityof(μs, 1) == 1.0
    @test densityof(μs, 0.5) == 0

    μs = SuperpositionMeasure([μ, ν])
    @test μs isa SuperpositionMeasure{<:AbstractVector{<:AbstractMeasure}}
    @test exp(logdensity_def(μs, 0)) == 1.0
    @test basemeasure(μs) == MeasureJoin([μ, ν])
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
    # Spike and slab with unnormalized weights:
    sm2 = 0.3 * N + 0.5 * Dirac(0.0)
    @test massof(sm2) ≈ 0.8
    @test logdensityof(sm2, 0.0) ≈ log(0.5)
    @test logdensityof(sm2, 1.0) ≈ log(0.3) + ld_N(1.0)
    @test logdensities(sm2, [0.0, 1.0]) ≈ [log(0.5), log(0.3) + ld_N(1.0)]
    draws = [rand(sm2) for _ in 1:10^5]
    @test isapprox(count(iszero, draws) / length(draws), 0.5 / 0.8, atol = 0.01)
    @test isapprox(count(iszero, rand(sm2^(10^5))) / 10^5, 0.5 / 0.8, atol = 0.01)

    @test MeasureBase.insupportof(MeasureBase.TrivialMeasure(), 1.0) == false
    @test logdensity_rel(MeasureBase.TrivialMeasure(), MeasureBase.LebesgueBase(), 1.0) == -Inf

    # Dirac equality is not decidable from types, so no weighted collapse,
    # the equal atoms join to one root:
    μ2 = μ + μ
    @test μ2 isa SuperpositionMeasure
    @test μ2 == superpose(μ, μ)
    @test exp(logdensity_def(μ2, 0)) == 2.0

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

@testset "roots of superpositions" begin
    pdfsum(ds, x) = log(sum(d -> pdf(d, x), ds))

    # Components with equal roots have the sum of their densities, in all
    # representations of the superposition:
    ds = (MvNormal([0.0, 1.0], [1.0 0.2; 0.2 1.5]), MvNormal([1.0, -1.0], [2.0 0.0; 0.0 0.5]))
    a, b = map(asmeasure, ds)
    x = [0.3, 0.4]
    for s in (superpose(a, b), a + b, superpose((a, b)), superpose([a, b]))
        @test logdensityof(s, x) ≈ pdfsum(ds, x)
        @test logdensity_rel(s, Lebesgue()^2, x) ≈ pdfsum(ds, x)
    end
    @test logdensityof(superpose(a, b, a), x) ≈ log(2 * pdf(ds[1], x) + pdf(ds[2], x))
    @test logdensityof(0.3 * a + 0.7 * b, x) ≈ log(0.3 * pdf(ds[1], x) + 0.7 * pdf(ds[2], x))
    ℓ_sum(x) = log(exp(logdensityof(StdNormal()^2, x)) + exp(logdensityof(StdExponential()^2, x)))
    @test logdensityof(StdNormal()^2 + StdExponential()^2, x) ≈ ℓ_sum(x)
    @test logdensityof(StdNormal()^static(2) + StdExponential()^static(2), SVector(0.3, 0.4)) ≈ ℓ_sum(x)
    dd = (Dirichlet([1.0, 2.0, 3.0]), Dirichlet([2.0, 2.0, 1.0]))
    @test logdensityof(asmeasure(dd[1]) + asmeasure(dd[2]), [0.2, 0.3, 0.5]) ≈ pdfsum(dd, [0.2, 0.3, 0.5])
    P2, P3 = asmeasure(Poisson(2.0))^2, asmeasure(Poisson(3.0))^2
    @test logdensityof(P2 + P3, [1, 2]) ≈ log(prod(pdf.(Poisson(2.0), [1, 2])) + prod(pdf.(Poisson(3.0), [1, 2])))

    # Mutually singular roots join to their sum, a point mass added to a
    # continuous measure has the weight of the atom at the atom:
    s = StdNormal() + 0.3 * Dirac(0.0)
    J = rootmeasure(s)
    @test J == MeasureJoin((LebesgueBase(), Dirac(0.0)))
    @test logdensityof(s, 0.0) ≈ log(0.3)
    @test logdensityof(s, 0.5) ≈ logpdf(Normal(), 0.5)
    @test logdensity_rel(J, Lebesgue(), 0.5) == 0 && logdensity_rel(J, Lebesgue(), 0.0) == Inf
    @test logdensity_rel(Lebesgue(), J, 0.0) == -Inf && logdensity_rel(Dirac(0.0), J, 0.0) == 0
    @test logdensity_rel(s, J, 0.0) ≈ log(0.3)

    # Joins relate to joins and to other superpositions:
    s2 = StdNormal() + Dirac(0.0) + Dirac(1.0)
    J2 = rootmeasure(s2)
    @test logdensity_rel(J, J2, 0.0) == 0 && logdensity_rel(J, J2, 1.0) == -Inf && logdensity_rel(J2, J, 1.0) == Inf
    @test logdensity_rel(J, J2, 0.5) == 0
    @test logdensity_rel(s, s2, 0.0) ≈ log(0.3)
    @test logdensity_rel(s, s2, 0.5) ≈ 0 atol = 1e-15
    @test logdensity_rel(s2, s, 1.0) == Inf

    # Roots of identical singleton types are the root itself:
    @test rootmeasure(StdNormal() + StdExponential()) === LebesgueBase()
    @test rootmeasure(superpose([StdNormal(), 2 * StdNormal()])) === LebesgueBase()
    @test logdensityof(superpose([StdNormal(), 2 * StdNormal()]), 0.3) ≈ log(3) + logpdf(Normal(), 0.3)

    # Joins of arrays of roots record whether all roots are equal:
    @test MeasureJoin([Dirac(0.0), Dirac(0.0)]).equal_roots
    @test !MeasureJoin([Dirac(0.0), Dirac(1.0)]).equal_roots
    let s = superpose([Dirac(0.0), 0.5 * Dirac(1.0), 0.25 * Dirac(0.0)])
        @test logdensityof(s, 0.0) ≈ log(1.25) && logdensityof(s, 1.0) ≈ log(0.5) && logdensityof(s, 0.5) == -Inf
    end
    @test logdensityof(superpose([Dirac(0.0), 2 * Dirac(0.0)]), 0.0) ≈ log(3)

    # Mixture models convert to superpositions with the pdf of the mixture:
    mix = Distributions.MixtureModel([MvNormal([0.0, 0.0], [1.0 0.0; 0.0 1.0]), MvNormal([1.0, 1.0], [2.0 0.0; 0.0 2.0])], [0.4, 0.6])
    @test logdensityof(asmeasure(mix), x) ≈ logpdf(mix, x)
end
