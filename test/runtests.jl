using Test
using Base.Iterators: take
using Random
using LinearAlgebra

using MeasureBase
using MeasureBase: logdensityof, logdensity_def

using Aqua
Aqua.test_all(MeasureBase; ambiguities = false, unbound_args = false)

d = ∫exp(x -> -x^2, Lebesgue(ℝ))

# if VERSION ≥ v"1.6"
#     @eval using JETTest

#     @eval begin     
#         @test_nodispatch density_def(Lebesgue(ℝ), 0.3)

#         @test_nodispatch density_def(Dirac(0), 0.3)
#         @test_nodispatch density_def(Dirac(0), 0)

#         @test_nodispatch density_def(d, 3)

#         @test_nodispatch basemeasure(d)

#         @test_nodispatch logdensity_def(For(3) do j Dirac(j) end, [1,2,3])
#     end
# end

# function draw2(μ)
#     x = rand(μ)
#     y = rand(μ)
#     while x == y
#         y = rand(μ)
#     end
#     return (x,y)
# end

function test_testvalue(μ)
    logdensity_def(μ, testvalue(μ)) isa AbstractFloat
end

test_measures = [
    # Chain(x -> Normal(μ=x), Normal(μ=0.0))
    For(3) do j
        Dirac(j)
    end
    For(2, 3) do i, j
        Dirac(i) + Dirac(j)
    end
    Lebesgue(ℝ)^3
    Lebesgue(ℝ)^(2, 3)
    3 * Lebesgue(ℝ)
    Dirac(π)
    Lebesgue(ℝ)
    Dirac(0.0) + Lebesgue(ℝ)
    SpikeMixture(Lebesgue(ℝ), 2)
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
        @test test_testvalue(μ)
    end

    for μ in testbroken_measures
        @info "testing $μ"
        @test_broken test_testvalue(μ)
    end
end

# @testset "Kernel" begin
#     κ = MeasureBase.kernel(MeasureBase.Dirac, identity)
#     @test rand(κ(1.1)) == 1.1
# end

@testset "SpikeMixture" begin
    @test rand(SpikeMixture(Dirac(0), 0.5)) == 0
    @test rand(SpikeMixture(Dirac(1), 1.0)) == 1
    w = 1 / 3
    m = SpikeMixture(d, w)
    bm = basemeasure(m)
    @test (bm.s * bm.w) * bm.m == 1.0 * basemeasure(d)
    @test density_def(m, 1.0) * (bm.s * bm.w) ≈ w * density_def(d, 1.0)
    @test density_def(m, 0) * (bm.s * (1 - bm.w)) ≈ (1 - w)
end

@testset "Dirac" begin
    @test rand(Dirac(0.2)) == 0.2
    @test logdensity_def(Dirac(0.3), 0.3) == 0.0
    @test logdensity_def(Dirac(0.3), 0.4) == -Inf
end

@testset "AffineTransform" begin
    f = AffineTransform((μ = 3, σ = 2))
    @test f(inv(f)(1)) == 1
    @test inv(f)(f(1)) == 1

    f = AffineTransform((μ = 3, ω = 2))
    @test f(inv(f)(1)) == 1
    @test inv(f)(f(1)) == 1

    f = AffineTransform((σ = 2,))
    @test f(inv(f)(1)) == 1
    @test inv(f)(f(1)) == 1

    f = AffineTransform((ω = 2,))
    @test f(inv(f)(1)) == 1
    @test inv(f)(f(1)) == 1

    f = AffineTransform((μ = 3,))
    @test f(inv(f)(1)) == 1
    @test inv(f)(f(1)) == 1
end

@testset "Affine" begin
    unif = ∫(x -> 0 < x < 1, Lebesgue(ℝ))
    f1 = AffineTransform((μ = 3.0, σ = 2.0))
    f2 = AffineTransform((μ = 3.0, ω = 2.0))
    f3 = AffineTransform((μ = 3.0,))
    f4 = AffineTransform((σ = 2.0,))
    f5 = AffineTransform((ω = 2.0,))

    for f in [f1, f2, f3, f4, f5]
        par = getfield(f, :par)
        @test Affine(par)(unif) == Affine(f, unif)
        @test density_def(Affine(f, Affine(inv(f), unif)), 0.5) == 1
    end

    d = ∫exp(x -> -x^2, Lebesgue(ℝ))

    μ = randn(3)
    σ = LowerTriangular(randn(3, 3))
    ω = inv(σ)

    x = randn(3)

    @test logdensity_def(Affine((μ = μ, σ = σ), d^3), x) ≈
          logdensity_def(Affine((μ = μ, ω = ω), d^3), x)
    @test logdensity_def(Affine((σ = σ,), d^3), x) ≈ logdensity_def(Affine((ω = ω,), d^3), x)
    @test logdensity_def(Affine((μ = μ,), d^3), x) ≈ logdensity_def(d^3, x - μ)

    d = ∫exp(x -> -x^2, Lebesgue(ℝ))
    a = Affine((σ = [1 0]',), d^1)
    x = randn(2)
    y = randn(1)
    @test logdensityof(a, x) ≈ logdensityof(d, inv(a.f)(x)[1])
    @test logdensityof(a, a.f(y)) ≈ logdensityof(d^1, y)

    b = Affine((ω = [1 0]'',), d^1)
    @test logdensityof(b, x) ≈ logdensityof(d, inv(b.f)(x)[1])
    @test logdensityof(b, b.f(y)) ≈ logdensityof(d^1, y)
end

@testset "For" begin
    FORDISTS = [
        For(1:10) do j
            Dirac(j)
        end
        For(4, 3) do i, j
            Dirac(i) ⊗ Dirac(j)
        end
        For(1:4, 1:4) do i, j
            Dirac(i) ⊗ Dirac(j)
        end
        For(eachrow(rand(4, 2))) do x
            Dirac(x[1]) ⊗ Dirac(x[2])
        end
        For(rand(4), rand(4)) do i, j
            Dirac(i) ⊗ Dirac(j)
        end
    ]

    for d in FORDISTS
        @info "testing $d"
        @test logdensity_def(d, rand(d)) isa Float64
    end
end

@testset "Half" begin
    Normal() = ∫exp(x -> -0.5x^2, Lebesgue(ℝ))
    @half Normal
    @test logdensity_def(HalfNormal(), Lebesgue(ℝ), -0.2) == -Inf
    @test logdensity_def(HalfNormal(), 0.2) == logdensity_def(Normal(), 0.2)
    @test density_def(HalfNormal(), Lebesgue(ℝ), 0.2) ≈ 2 * density_def(Normal(), Lebesgue(ℝ), 0.2)
end

# import MeasureBase.:⋅
# function ⋅(μ::Normal, kernel) 
#     m = kernel(μ)
#     Normal(μ = m.μ.μ, σ = sqrt(m.μ.σ^2 + m.σ^2))
# end
# struct AffineMap{S,T}
#     B::S
#     β::T
# end
# (a::AffineMap)(x) = a.B*x + a.β
# (a::AffineMap)(p::Normal) = Normal(μ = a.B*mean(p) + a.β, σ = sqrt(a.B*p.σ^2*a.B'))

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
