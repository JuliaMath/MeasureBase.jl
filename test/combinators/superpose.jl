using Test

using MeasureBase
using MeasureBase: superpose, weightedmeasure, StdNormal

@testset "superpose.jl" begin
    μ = Dirac(0)
    ν = Dirac(1)
    μs = μ + ν
    @test μs isa SuperpositionMeasure{<:Tuple{Dirac,Dirac}}
    @test μs == SuperpositionMeasure((μ, ν)) == superpose(μ, ν)
    @test density_def(μs, 0) == 1.0
    @test basemeasure(μs) == CountingBase() + CountingBase()

    μs = SuperpositionMeasure([μ, ν])
    @test μs isa SuperpositionMeasure{<:AbstractVector{<:AbstractMeasure}}
    @test density_def(μs, 0) == 1.0
    @test basemeasure(μs) == weightedmeasure(log(2), CountingBase())

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
