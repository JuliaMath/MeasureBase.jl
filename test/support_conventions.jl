# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, StdLogistic, Half, logdensities
using MeasureBase: weightedmeasure, productmeasure, mcombine
using JLArrays: JLArray

# Variates of the right shape never throw: densities are -Inf and
# transports are NaN outside the support, wrong shapes throw.
@testset "support conventions" begin
    outside = Dict(StdUniform() => -0.5, StdExponential() => -1.0, Half(StdNormal()) => -1.0)

    @testset "densities" begin
        for (μ, x) in outside
            @test logdensityof(μ, x) == -Inf
            @test logdensityof(weightedmeasure(0.3, μ), x) == -Inf
            X = [x 0.5 0.5; 0.5 x 0.5; 0.5 0.5 0.5]
            ℓ = logdensities(μ^3, X)
            @test ℓ[1:2] == [-Inf, -Inf] && isfinite(ℓ[3])
            @test Array(logdensities(μ^3, JLArray(X))) == ℓ
            @test Array(logdensities(μ, JLArray(vec(X)))) == logdensities(μ, vec(X))
        end
        @test logdensityof(StdUniform(), Inf) == -Inf
        @test logdensityof(StdExponential(), Inf) == -Inf
        @test logdensityof(StdNormal(), -Inf) == -Inf
        @test logdensityof(productmeasure((StdUniform(), StdExponential())), (0.5, -1.0)) == -Inf
        @test logdensityof(mcombine(vcat, StdUniform()^2, StdExponential()^1), [0.5, 1.5, 0.5]) == -Inf
    end

    @testset "transports" begin
        stds = (StdUniform(), StdExponential(), StdLogistic(), StdNormal())
        for (μ, x) in outside, ν in stds
            μ === ν && continue
            f = transport_to(ν, μ)
            @test isnan(f(x))
            @test isnan(transport_to(ν^2, μ^2)([x, 0.5])[1])
            Y = f.([x, 0.5, 0.5])
            @test isnan(Y[1]) && !isnan(Y[2])
            @test isequal(Array(f.(JLArray([x, 0.5, 0.5]))), Y)
            @test isnan(transport_to(μ, ν)(NaN))
        end
        for ν in (StdExponential(), StdLogistic(), StdNormal(), Half(StdNormal()))
            @test isnan(transport_to(ν, StdUniform())(1.5))
            @test isnan(transport_to(ν, StdUniform())(-0.5))
        end
        # Endpoints of the unit interval stand for their nearest interior
        # points, tails never underflow to infinite variates:
        for ν in (StdExponential(), StdLogistic(), StdNormal(), Half(StdNormal()))
            f = transport_to(ν, StdUniform())
            @test isfinite(f(0.0)) && isfinite(f(1.0)) && f(0.0) <= f(0.5) <= f(1.0)
            @test f(1.0) == f(prevfloat(1.0)) && f(0.0) == f(floatmin(Float64))
        end
        for (ν, μ) in ((StdExponential(), StdNormal()), (StdNormal(), StdExponential()), (StdNormal(), StdLogistic()))
            @test all(isfinite, transport_to(ν, μ).([-1e6, -40.0, 40.0, 1e6][MeasureBase.insupport.(Ref(μ), [-1e6, -40.0, 40.0, 1e6])]))
        end
        @test !isnan(transport_to(StdUniform(), StdNormal())(-37.0))
        @test !isnan(transport_to(StdUniform(), StdLogistic())(-800.0))
    end

    @testset "wrong shapes throw" begin
        @test_throws ArgumentError logdensityof(StdNormal()^3, randn(2))
        @test_throws ArgumentError logdensities(StdNormal()^3, randn(2, 5))
        @test_throws ArgumentError transport_to(StdUniform()^3, StdNormal()^3)(randn(2))
        @test_throws ArgumentError logdensityof(StdNormal(), randn(2))
    end
end
