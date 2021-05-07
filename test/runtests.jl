using MeasureBase
using Test
using Base.Iterators: take
using Random
using LinearAlgebra
using KeywordCalls
using Statistics


function draw2(μ)
    x = rand(μ)
    y = rand(μ)
    while x == y
        y = rand(μ)
    end
    return (x,y)
end

const sqrt2π = sqrt(2π)

@testset "Parameterized Measures" begin
    @measure Normal(μ,σ)
    @kwstruct Normal(μ)
    @kwstruct Normal()
    
    MeasureBase.basemeasure(::Normal)= (1/sqrt2π) * Lebesgue(ℝ)
    MeasureBase.logdensity(d::Normal{(:μ,:σ)}, x) = -log(d.σ) - (x - d.μ)^2 / (2 * d.σ^2)
    MeasureBase.logdensity(d::Normal{(:μ,)}, x) = - (x - d.μ)^2 / 2
    MeasureBase.logdensity(d::Normal{()}, x) = - x^2 / 2

    Base.rand(rng::Random.AbstractRNG, T::Type, d::Normal{(:μ,:σ)}) = d.μ + d.σ * randn(rng, T)
    Base.rand(rng::Random.AbstractRNG, T::Type, d::Normal{(:μ,)}) = d.μ + randn(rng, T)
    Base.rand(rng::Random.AbstractRNG, T::Type, d::Normal{()}) = randn(rng, T)

    MeasureBase.representative(d::Normal{(:μ,:σ)}) = d.σ > 0.0 ? Lebesgue(ℝ) : Dirac(d.μ)
    MeasureBase.representative(d::Normal{(:μ,)}) = Lebesgue(ℝ)

    # Leave this undefined to test fallback inference algorithm
    # MeasureBase.representative(::Normal) = Lebesgue(ℝ)

    @test Normal(2,4) == Normal(μ=2, σ=4)
    @test Normal(σ=4, μ=2) == Normal(μ=2, σ=4)
    @test logdensity(Normal(), 3) == logdensity(Normal(0,1), 3)

    x = randn()
    @test_broken logdensity(Normal(3,2), Lebesgue(ℝ), x) ≈ logdensity(Normal(3,2), Normal(), x ) + logdensity(Normal(), Lebesgue(ℝ),x)
    @test_broken 𝒹(Normal(3,2), Normal())(x) ≈ logdensity(Normal(3,2), Normal(), x)
end

@testset "Density" begin
    x = randn()
    f(x) = -x^2
    μ = Normal()
    ν = Lebesgue(ℝ)
    @test_broken 𝒹(∫(f, μ), μ)(x) ≈ f(x)
    @test_broken logdensity(∫(𝒹(μ, ν), ν), x) ≈ logdensity(μ, x)
end


@testset "Kernel" begin
    κ = kernel(identity, Dirac)
    @test rand(κ(1.1)) == 1.1
    @test kernelize(Normal(0,1)) == (Kernel{Normal, UnionAll}(NamedTuple{(:μ, :σ), T} where T<:Tuple), (0, 1))
end

@testset "SpikeMixture" begin
    @test rand(SpikeMixture(Dirac(0), 0.5)) == 0
    @test rand(SpikeMixture(Dirac(1), 1.0)) == 1
    w = 1/3
    m = SpikeMixture(Normal(), w)
    bm = basemeasure(m)
    @test (bm.s*bm.w)*bm.m == 1.0*basemeasure(Normal())
    @test density(m, 1.0)*(bm.s*bm.w) == w*density(Normal(),1.0)
    @test density(m, 0)*(bm.s*(1-bm.w)) ≈ (1-w)
end

@testset "For" begin
    FORDISTS = [
        For(1:10) do j Normal(μ=j) end
        For(4,3) do μ,σ Normal(μ,σ) end
        For(1:4, 1:4) do μ,σ Normal(μ,σ) end
        For(eachrow(rand(4,2))) do x Normal(x[1], x[2]) end
        For(rand(4), rand(4)) do μ,σ Normal(μ,σ) end
    ]

    for d in FORDISTS
        @test logdensity(d, rand(d)) isa Float64
    end
end

function ⋅(μ::Normal, kernel) 
    m = kernel(μ)
    Normal(μ = m.μ.μ, σ = sqrt(m.μ.σ^2 + m.σ^2))
end

"""
    ConstantMap(β)
Represents a function `f = ConstantMap(β)`
such that `f(x) == β`.
"""
struct ConstantMap{T}
    x::T
end
(a::ConstantMap)(x) = a.x
(a::ConstantMap)() = a.x

struct AffineMap{S,T}
    B::S
    β::T
end
(a::AffineMap)(x) = a.B*x + a.β
(a::AffineMap)(p::Normal) = Normal(μ = a.B*mean(p) + a.β, σ = sqrt(a.B*p.σ^2*a.B'))

@testset "DynamicFor" begin
    mc = Chain(Normal(μ=0.0)) do x Normal(μ=x) end
    r = rand(mc)
   
    # Check that `r` is now deterministic
    @test logdensity(mc, take(r, 100)) == logdensity(mc, take(r, 100))
    
    d2 = For(r) do x Normal(μ=x) end  

    @test_broken let r2 = rand(d2)
        logdensity(d2, take(r2, 100)) == logdensity(d2, take(r2, 100))
    end
end

@testset "LogLikelihood" begin
    d = Normal()
    ℓ = LogLikelihood(Normal{(:μ,)}, 3.0) 
    @test logdensity(d ⊙ ℓ, 2.0) == logdensity(d, 2.0) + logdensity(ℓ, 2.0)
end
