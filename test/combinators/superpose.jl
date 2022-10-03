using Test

using MeasureBase
using MeasureBase: superpose

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
    @test_throws ErrorException density_def(μs, 0)
    @test basemeasure(μs).components ==
          SuperpositionMeasure([CountingBase(), CountingBase()]).components

    μ2 = μ + μ
    @test μ2 isa WeightedMeasure
    @test μ2 == superpose(μ, μ)
    @test basemeasure(μ2) == μ
end
