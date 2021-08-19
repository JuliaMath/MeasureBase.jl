using Test
using Base.Iterators: take
using Random
using LinearAlgebra

using MeasureBase

using Aqua
Aqua.test_all(MeasureBase; ambiguities=false, unbound_args=false)

# function draw2(μ)
#     x = rand(μ)
#     y = rand(μ)
#     while x == y
#         y = rand(μ)
#     end
#     return (x,y)
# end

# function test_measure(μ)
#     logdensity(μ, testvalue(μ)) isa AbstractFloat
# end

# test_measures = [
#     # Chain(x -> Normal(μ=x), Normal(μ=0.0))
#     For(3) do j Normal(σ=j) end
#     For(2,3) do i,j Normal(i,j) end
#     Lebesgue(ℝ) ^ 3
#     Lebesgue(ℝ) ^ (2,3)
#     3 * Lebesgue(ℝ)
#     Dirac(π)
#     Lebesgue(ℝ)
#     # Normal() ⊙ Cauchy()
# ]

# testbroken_measures = [
#     Pushforward(as𝕀, Normal())
#     SpikeMixture(Normal(), 2)
#     # InverseGamma(2) # Not defined yet
#     # MvNormal(I(3)) # Entirely broken for now
#     CountingMeasure(Float64)
#     Likelihood
#     Dirac(0.0) + Lebesgue(ℝ)

#     TrivialMeasure()
# ]

# @testset "testvalue" begin
#     for μ in test_measures
#         @test test_measure(μ)
#     end

#     for μ in testbroken_measures
#         @test_broken test_measure(μ)
#     end
    
#     @testset "testvalue(::Chain)" begin
#         mc =  Chain(x -> Normal(μ=x), Normal(μ=0.0))
#         r = testvalue(mc)
#         @test logdensity(mc, Iterators.take(r, 10)) isa AbstractFloat
#     end
# end


# @testset "Kernel" begin
#     κ = MeasureBase.kernel(MeasureBase.Dirac, identity)
#     @test rand(κ(1.1)) == 1.1
# end

# @testset "SpikeMixture" begin
#     @test rand(SpikeMixture(Dirac(0), 0.5)) == 0
#     @test rand(SpikeMixture(Dirac(1), 1.0)) == 1
#     w = 1/3
#     m = SpikeMixture(Normal(), w)
#     bm = basemeasure(m)
#     @test (bm.s*bm.w)*bm.m == 1.0*basemeasure(Normal())
#     @test density(m, 1.0)*(bm.s*bm.w) == w*density(Normal(),1.0)
#     @test density(m, 0)*(bm.s*(1-bm.w)) ≈ (1-w)
# end

# @testset "Dirac" begin
#     @test rand(Dirac(0.2)) == 0.2
#     @test logdensity(Dirac(0.3), 0.3) == 0.0
#     @test logdensity(Dirac(0.3), 0.4) == -Inf
# end

# @testset "For" begin
#     FORDISTS = [
#         For(1:10) do j Normal(μ=j) end
#         For(4,3) do μ,σ Normal(μ,σ) end
#         For(1:4, 1:4) do μ,σ Normal(μ,σ) end
#         For(eachrow(rand(4,2))) do x Normal(x[1], x[2]) end
#         For(rand(4), rand(4)) do μ,σ Normal(μ,σ) end
#     ]

#     for d in FORDISTS
#         @test logdensity(d, rand(d)) isa Float64
#     end
# end

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

# @testset "DynamicFor" begin
#     mc = Chain(Normal(μ=0.0)) do x Normal(μ=x) end
#     r = rand(mc)
   
#     # Check that `r` is now deterministic
#     @test logdensity(mc, take(r, 100)) == logdensity(mc, take(r, 100))
    
#     d2 = For(r) do x Normal(μ=x) end  

#     @test_broken let r2 = rand(d2)
#         logdensity(d2, take(r2, 100)) == logdensity(d2, take(r2, 100))
#     end
# end


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
#             @test logdensity(d ⊙ ℓ, p) == logdensity(d, p) + logdensity(ℓ, p)
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

# @testset "Show methods" begin
#     @testset "PowerMeasure" begin
#         @test repr(Lebesgue(ℝ) ^ 5) == "Lebesgue(ℝ) ^ 5"
#         @test repr(Lebesgue(ℝ) ^ (3, 2)) == "Lebesgue(ℝ) ^ (3, 2)"
#     end
# end

# @testset "Density measures and Radon-Nikodym" begin
#     x = randn()
#     let d = ∫(𝒹(Cauchy(), Normal()), Normal())
#         @test logdensity(d, x) ≈ logdensity(Cauchy(), x) 
#     end

#     let f = 𝒹(∫(x -> x^2, Normal()), Normal())
#         @test f(x) ≈ x^2
#     end

#     let d = ∫exp(log𝒹(Cauchy(), Normal()), Normal())
#         @test logdensity(d, x) ≈ logdensity(Cauchy(), x) 
#     end

#     let f = log𝒹(∫exp(x -> x^2, Normal()), Normal())
#         @test f(x) ≈ x^2
#     end
# end
