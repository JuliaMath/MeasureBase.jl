# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, logdensities, weightedmeasure, productmeasure, pushfwd
using MeasureBase: batched_transport_to_std, batched_transport_from_std
using InverseFunctions: inverse
using FixedSizeArrays: FixedSizeArrayDefault
using ArraysOfArrays: flatview, sliced
using AffineMaps: Mul
using Distributions: Normal

# Fixed-size inputs give fixed-size outputs, the kernels allocate via
# `similar` and never fall back to plain arrays:
@testset "fixed-size arrays" begin
    fixed(A) = FixedSizeArrayDefault(A)
    isfixed(A) = A isa FixedSizeArrayDefault
    isfixed(A::Union{SubArray,Base.ReshapedArray}) = isfixed(parent(A))
    X = fixed(randn(3, 20))
    x = fixed(randn(3))
    m3 = StdNormal()^3

    @testset "densities" begin
        ℓ = logdensities(m3, X)
        @test isfixed(ℓ) && ℓ ≈ logdensities(m3, Array(X))
        @test logdensityof(m3, x) ≈ logdensityof(m3, Array(x))
        @test isfixed(logdensities(m3, sliced(X, Val(1))))
        @test isfixed(logdensities((StdNormal()^2)^3, fixed(randn(2, 3, 4))))
        @test isfixed(logdensities(weightedmeasure(0.3, m3), X))
        P = productmeasure(fixed([pushfwd(Mul(s), StdNormal()) for s in (1.0, 2.0, 3.0)]))
        @test isfixed(logdensities(P, X))
        Pn = productmeasure(fixed([Normal(μ, 1.0) for μ in (0.0, 1.0, 2.0)]))
        @test isfixed(logdensities(Pn, X)) && logdensities(Pn, X) ≈ logdensities(Pn, Array(X))
    end

    @testset "transports" begin
        f = transport_to(StdUniform()^3, m3)
        y = f(x)
        @test isfixed(y) && y ≈ f(Array(x))
        Y = f.(X)
        @test isfixed(flatview(Y)) && flatview(Y) ≈ flatview(f.(Array(X)))
        @test flatview(inverse(f).(Y)) ≈ X
        @test isfixed(batched_transport_to_std(StdNormal, StdExponential()^3, fixed(rand(3, 5))))
        @test isfixed(batched_transport_from_std(StdNormal, StdExponential()^3, fixed(randn(3, 5))))
        @test isfixed(MeasureBase.convert_realtype(Float32, x))
    end
end
