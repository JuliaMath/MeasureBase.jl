# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Densities, derivatives and tails of the wrapped distribution families
# against independent references: analytic values and higher precision.

using Test

using LinearAlgebra
using Distributions
import AffineMaps
import ForwardDiff, Zygote, Mooncake

using MeasureBase
using MeasureBase: StdNormal, StdUniform, transport_to, logdensities

@testset "reference values" begin
    @testset "densities of affine pushforwards" begin
        # The law of `2z` for a standard normal `z`:
        p = pushfwd(AffineMaps.Mul(2.0I(2)), StdNormal()^2)
        Y = [0.0 0.3 -1.2; 0.0 1.1 0.4]
        ℓ_ref = logpdf(MvNormal(zeros(2), 4.0I(2)), Y)
        @test logdensity_def(p, Y[:, 2]) ≈ ℓ_ref[2]
        @test MeasureBase.batched_logdensity_def(p, zeros(2)) ≈ ℓ_ref[1]
        @test collect(MeasureBase.batched_logdensity_def(p, Y)) ≈ ℓ_ref
        @test logdensities(p, Y) ≈ ℓ_ref
        @test logdensity_def(p^3, Y) ≈ sum(ℓ_ref)
        @test logdensity_rel(p^3, basemeasure(p)^3, Y) ≈ sum(ℓ_ref)
        @test logdensityof(p^3, Y) ≈ sum(ℓ_ref)
    end
end
