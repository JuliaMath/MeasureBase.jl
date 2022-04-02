using Test
using Base.Iterators: take
using Random
using LinearAlgebra

using MeasureBase
using MeasureBase: test_interface

using Aqua
Aqua.test_all(MeasureBase; ambiguities = false, unbound_args = false)

d = ∫exp(x -> -x^2, Lebesgue(ℝ))

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
    # Dirac(0.0) + Lebesgue(ℝ)
    SpikeMixture(Lebesgue(ℝ), 0.2)
    # d ⊙ d
]

testbroken_measures = [
    # InverseGamma(2) # Not defined yet
    # MvNormal(I(3)) # Entirely broken for now
    CountingMeasure()
    Likelihood
    TrivialMeasure()
]

@testset "testvalue" begin
    for μ in test_measures
        @info "testing $μ"
        test_interface(μ)
        test_interface(μ ^ 3)
        test_interface(μ ^ (3,2))
        test_interface(5 * μ)
        # test_interface(SpikeMixture(μ, 0.2))
    end

    for μ in testbroken_measures
        @info "testing $μ"
        @test_broken test_testvalue(μ)
    end
end

# @testset "Kleisli" begin
#     κ = MeasureBase.kleisli(MeasureBase.Dirac, identity)
#     @test rand(κ(1.1)) == 1.1
# end

@testset "SpikeMixture" begin
    @test rand(SpikeMixture(Dirac(0), 0.5)) == 0
    @test rand(SpikeMixture(Dirac(1), 1.0)) == 1
end

@testset "Dirac" begin
    @test rand(Dirac(0.2)) == 0.2
    @test logdensity_def(Dirac(0.3), 0.3) == 0.0
    @test logdensity_def(Dirac(0.3), 0.4) == -Inf
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

@testset "Half" begin
    Normal() = ∫exp(x -> -0.5x^2, Lebesgue(ℝ))
    HalfNormal() = Half(Normal())
    @test logdensityof(HalfNormal(), -0.2) == -Inf
    @test logdensity_def(HalfNormal(), 0.2) == logdensity_def(Normal(), 0.2)
    @test densityof(HalfNormal(), 0.2) ≈ 2 * densityof(Normal(), 0.2)
end

@testset "Likelihood" begin
    ℓ = Likelihood(3) do (μ,)
        ∫exp(Lebesgue(ℝ)) do x
            -(x - μ)^2
        end
    end

    @inferred logdensityof(ℓ, 2.0)
end

# @testset "Likelihood" begin
#     dps = [
#         (Normal()                             ,    2.0  )
#         # (Pushforward(as((μ=asℝ,)), Normal()^1), (μ=2.0,))
#     ]

#     ℓs = [
#         Likelihood(Normal{(:μ,)},              3.0)
#         Likelihood(kleisli(Normal, x -> (μ=x, σ=2.0)), 3.0)
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
    @test MeasureBase.logdensity_rel(Dirac(0.0)+Lebesgue(), Dirac(0.0), 0.0) == 0.0
    @test MeasureBase.logdensity_rel(Dirac(0.0)+Lebesgue(), Dirac(0.0), 1.0) == Inf
    @test MeasureBase.logdensity_rel(Dirac(0.0)+Lebesgue(), Lebesgue(), 0.0) == Inf
    @test MeasureBase.logdensity_rel(Dirac(0.0)+Lebesgue(), Lebesgue(), 1.0) == 0.0
    @test MeasureBase.logdensity_rel(Lebesgue()+ Dirac(0.0), Dirac(0.0), 0.0) == 0.0
    @test MeasureBase.logdensity_rel(Lebesgue()+ Dirac(0.0), Dirac(0.0), 1.0) == Inf
    @test MeasureBase.logdensity_rel(Lebesgue()+ Dirac(0.0), Lebesgue(), 0.0) == Inf
    @test MeasureBase.logdensity_rel(Lebesgue()+ Dirac(0.0), Lebesgue(), 1.0) == 0.0
end

# @testset "Density measures and Radon-Nikodym" begin
#     x = randn()
#     let d = ∫(𝒹(Cauchy(), Normal()), Normal())
#         @test logdensity_def(d, x) ≈ logdensity_def(Cauchy(), x) 
#     end

#     let f = 𝒹(∫(x -> x^2, Normal()), Normal())
#         @test f(x) ≈ x^2
#     end

#     let d = ∫exp(log𝒹(Cauchy(), Normal()), Normal())
#         @test logdensity_def(d, x) ≈ logdensity_def(Cauchy(), x) 
#     end

#     let f = log𝒹(∫exp(x -> x^2, Normal()), Normal())
#         @test f(x) ≈ x^2
#     end
# end
