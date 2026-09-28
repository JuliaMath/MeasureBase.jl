using Test

using MeasureBase
using MeasureBase: superpose, weightedmeasure, StdNormal, StdUniform, SpikeMixture

@testset "superpose.jl" begin
    μ = Dirac(0)
    ν = Dirac(1)
    μs = μ + ν
    @test μs isa SuperpositionMeasure{<:Tuple{Dirac,Dirac}}
    @test μs == SuperpositionMeasure((μ, ν)) == superpose(μ, ν)
    @test density_def(μs, 0) == 1.0
    @test basemeasure(μs) == μs
    @test densityof(μs, 0) == 1.0
    @test densityof(μs, 1) == 1.0
    @test densityof(μs, 0.5) == 0

    μs = SuperpositionMeasure([μ, ν])
    @test μs isa SuperpositionMeasure{<:AbstractVector{<:AbstractMeasure}}
    @test density_def(μs, 0) == 1.0
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

    # Dirac equality is not decidable from types, so no weighted collapse:
    μ2 = μ + μ
    @test μ2 isa SuperpositionMeasure
    @test μ2 == superpose(μ, μ)
    @test density_def(μ2, 0) == 1.0

    # For singleton measure types equal measures combine into weighted measures:
    s2 = StdNormal() + StdNormal()
    @test s2 isa WeightedMeasure
    @test exp(s2.logweight) ≈ 2
    @test basemeasure(s2) == StdNormal()
end
