using Test
using Base.Iterators: take
using Random
using LinearAlgebra
import LogarithmicNumbers

using MeasureBase
using MeasureBase: test_interface, test_smf

using Aqua
Aqua.test_all(MeasureBase; ambiguities = false)

include("static.jl")

# Aqua._test_ambiguities(
#     Aqua.aspkgids(MeasureBase);
#     exclude = [LogarithmicNumbers.Logarithmic],
#     # packages::Vector{PkgId};
#     # color::Union{Bool, Nothing} = nothing,
#     # exclude::AbstractArray = [],
#     # # Options to be passed to `Test.detect_ambiguities`:
#     # detect_ambiguities_options...,
# )

d = mintegrate_exp(x -> -x^2, Lebesgue(ℝ))

# function draw2(μ)
#     x = rand(μ)
#     y = rand(μ)
#     while x == y
#         y = rand(μ)
#     end
#     return (x,y)
# end

test_measures = [
    # Chain(x -> Normal(μ=x), Normal(μ=0.0))
    # For(3) do j
    #     Dirac(j)
    # end
    # For(2, 3) do i, j
    #     Dirac(i) + Dirac(j)
    # end
    Lebesgue(ℝ)^3
    Lebesgue(ℝ)^(2, 3)
    3 * Lebesgue(ℝ)
    Dirac(π)
    Lebesgue(ℝ)
    0.2 * Lebesgue(ℝ) + 0.8 * Dirac(0.0)
    Dirac(0) + Dirac(1)
    Dirac(0.0) + Lebesgue(ℝ)
    SpikeMixture(Lebesgue(ℝ), 0.2)
    StdLogistic()
    StdLogistic()^3
    StdLogistic()^(2, 3)
    3 * StdLogistic()
    0.2 * StdLogistic() + 0.8 * Dirac(0.0)
    Dirac(0.0) + StdLogistic()
    SpikeMixture(StdLogistic(), 0.2)
    StdUniform()
    StdUniform()^3
    StdUniform()^(2, 3)
    3 * StdUniform()
    0.2 * StdUniform() + 0.8 * Dirac(0.0)
    Dirac(0.0) + StdUniform()
    SpikeMixture(StdUniform(), 0.2)
    StdExponential()^3
    StdExponential()^(2, 3)
    3 * StdExponential()
    StdExponential()
    0.2 * StdExponential() + 0.8 * Dirac(0.0)
    Dirac(0.0) + StdExponential()
    SpikeMixture(StdExponential(), 0.2)

    # d ⊙ d
]

testbroken_measures = [
    # InverseGamma(2) # Not defined yet
    # MvNormal(I(3)) # Entirely broken for now
    CountingBase()
    Likelihood
    TrivialMeasure()
]

@testset "testvalue" begin
    for μ in test_measures
        @info "testing $μ"
        test_interface(μ)
        test_interface(μ^3)
        test_interface(μ^(3, 2))
        test_interface(5 * μ)
        # test_interface(SpikeMixture(μ, 0.2))
    end

    for μ in testbroken_measures
        @info "testing $μ"
        @test_broken test_testvalue(μ)
    end
end

# @testset "TransitionKernel" begin
#     κ = MeasureBase.kernel(MeasureBase.Dirac, identity)
#     @test rand(κ(1.1)) == 1.1
# end

@testset "SpikeMixture" begin
    @test rand(SpikeMixture(Dirac(0), 0.5)) == 0
    @test rand(SpikeMixture(Dirac(1), 1.0)) == 1
end

@testset "Dirac" begin
    @test rand(Dirac(0.2)) == 0.2
    @test logdensityof(Dirac(0.3), 0.3) == 0.0
    @test logdensityof(Dirac(0.3), 0.4) == -Inf
end

# @testset "For" begin
#     FORDISTS = [
#         For(1:10) do j
#             Dirac(j)
#         end
#         For(4, 3) do i, j
#             Dirac(i) ⊗ Dirac(j)
#         end
#         For(1:4, 1:4) do i, j
#             Dirac(i) ⊗ Dirac(j)
#         end
#         For(eachrow(rand(4, 2))) do x
#             Dirac(x[1]) ⊗ Dirac(x[2])
#         end
#         For(rand(4), rand(4)) do i, j
#             Dirac(i) ⊗ Dirac(j)
#         end
#     ]

#     for d in FORDISTS
#         test_interface(d)
#     end
# end

@testset "broadcasting" begin
    @test logdensityof.(Dirac(2), [1, 2, 3]) isa Vector{Float64}
end

@testset "powers" begin
    @test @inferred productmeasure(fill(Lebesgue(), 5)) isa PowerMeasure
    @test logdensityof(Lebesgue()^3, 2) == logdensityof(Lebesgue()^(3,), 2)
    @test logdensityof(Lebesgue()^3, 2) == logdensityof(Lebesgue()^(3, 1), (2, 0))
end

Normal() = mintegrate_exp(x -> -0.5x^2, Lebesgue(ℝ))

@testset "Half" begin
    HalfNormal() = Half(Normal())
    @test logdensityof(HalfNormal(), -0.2) == -Inf
    @test logdensity_def(HalfNormal(), 0.2) == logdensity_def(Normal(), 0.2)
    @test densityof(HalfNormal(), 0.2) ≈ 2 * densityof(Normal(), 0.2)
end

@testset "Likelihood" begin
    ℓ = Likelihood(3) do (μ,)
        mintegrate_exp(Lebesgue(ℝ)) do x
            -(x - μ)^2
        end
    end
end

# @testset "Likelihood" begin
#     dps = [
#         (Normal()                             ,    2.0  )
#         # (Pushforward(as((μ=asℝ,)), Normal()^1), (μ=2.0,))
#     ]

#     ℓs = [
#         Likelihood(Normal{(:μ,)},              3.0)
#         Likelihood(kernel(Normal, x -> (μ=x, σ=2.0)), 3.0)
#     ]

#     for (d,p) in dps
#         for ℓ in ℓs
#             @test logdensity_def(d ⊙ ℓ, p) == logdensity_def(d, p) + logdensity_def(ℓ, p)
#         end
#     end
# end

# @testset "ProductMeasure" begin
#     d = For(1:10) do j Poisson(exp(j)) end
#     x = Vector{Int16}(undef, 10)
#     @test rand!(d,x) isa Vector
#     @test rand(d) isa Vector

#     @testset "Indexed by Generator" begin
#         d = For((j^2 for j in 1:10)) do i Poisson(i) end
#         x = Vector{Int16}(undef, 10)
#         @test rand!(d,x) isa Vector
#         @test_broken rand(d) isa Base.Generator
#     end

#     @testset "Indexed by multiple Ints" begin
#         d = For(2,3) do μ,σ Normal(μ,σ) end
#         x = Matrix{Float16}(undef, 2, 3)
#         @test rand!(d, x) isa Matrix
#         @test_broken rand(d) isa Matrix{Float16}
#     end
# end

@testset "Show methods" begin
    @testset "PowerMeasure" begin
        @test repr(Lebesgue(ℝ)^5) == "Lebesgue(ℝ) ^ 5"
        @test repr(Lebesgue(ℝ)^(3, 2)) == "Lebesgue(ℝ) ^ (3, 2)"
    end
end

@testset "logdensity_rel" begin
    @test logdensity_rel(Dirac(0.0) + Lebesgue(), Dirac(1.0), 0.0) == Inf
    @test logdensity_rel(Dirac(0.0) + Lebesgue(), Dirac(1.0), 1.0) == -Inf
    @test logdensity_rel(Dirac(0.0) + Lebesgue(), Dirac(1.0), 2.0) == Inf
    @test logdensity_rel(Dirac(0.0) + Lebesgue(), Dirac(0.0), 0.0) == 0.0
    @test logdensity_rel(Dirac(0.0) + Lebesgue(), Dirac(0.0), 1.0) == Inf
    @test logdensity_rel(Dirac(0.0) + Lebesgue(), Lebesgue(), 0.0) == Inf
    @test logdensity_rel(Dirac(0.0) + Lebesgue(), Lebesgue(), 1.0) == 0.0

    @test logdensity_rel(Dirac(1.0), Dirac(0.0) + Lebesgue(), 0.0) == -Inf
    @test logdensity_rel(Dirac(1.0), Dirac(0.0) + Lebesgue(), 1.0) == Inf
    @test logdensity_rel(Dirac(1.0), Dirac(0.0) + Lebesgue(), 2.0) == -Inf
    @test logdensity_rel(Dirac(0.0), Dirac(0.0) + Lebesgue(), 0.0) == 0.0
    @test logdensity_rel(Dirac(0.0), Dirac(0.0) + Lebesgue(), 1.0) == -Inf
    @test logdensity_rel(Lebesgue(), Dirac(0.0) + Lebesgue(), 0.0) == -Inf
    @test logdensity_rel(Lebesgue(), Dirac(0.0) + Lebesgue(), 1.0) == 0.0

    @test isnan(logdensity_rel(Dirac(0), Dirac(1), 2))
end

@testset "Density measures and Radon-Nikodym" begin
    x = randn()
    f(x) = x^2
    @test log(density_rel(mintegrate_exp(f, Lebesgue()), Lebesgue())(x)) ≈ f(x)

    let f = density_rel(mintegrate_exp(x -> x^2, Lebesgue()), Lebesgue())
        @test log(f(x)) ≈ x^2
    end

    let f = logdensity_rel(mintegrate_exp(x -> x^2, Normal()), Normal())
        @test f(x) ≈ x^2
    end
end

include("getdof.jl")
include("transport.jl")
include("smf.jl")

include("combinators/weighted.jl")
include("combinators/transformedmeasure.jl")
include("combinators/combined.jl")
include("combinators/bind.jl")

include("measure_operators.jl")
