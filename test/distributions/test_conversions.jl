# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random, Statistics, LinearAlgebra
using Distributions
using StableRNGs

import MeasureBase
using MeasureBase: AbstractMeasure, AsMeasure, asmeasure
using MeasureBase: StdUniform, StdNormal, StdExponential, StdLogistic
using MeasureBase: SuperpositionMeasure, PushforwardMeasure, ProductMeasure
using MeasureBase: logdensityof, massof, insupport


@testset "conversions" begin
    stblrng() = StableRNG(789990641)

    function test_conversion(d::Distribution, ::Type{M}) where {M}
        @testset "conversion $(typeof(d).name) <-> $M" begin
            m = asmeasure(d)
            @test m isa M
            @test typeof(convert(AbstractMeasure, d)) === typeof(m)
            @test_throws ArgumentError AsMeasure{typeof(d)}(d)

            d2 = convert(Distribution, m)
            @test d2 isa Distribution
            @test typeof(Distributions.Distribution(m)) === typeof(d2)

            for x in (rand(stblrng(), d) for _ in 1:10)
                @test logdensityof(m, x) ≈ logpdf(d, x)
                @test logpdf(d2, x) ≈ logpdf(d, x)
                @test insupport(m, x) != false
            end

            x = rand(stblrng(), Float64, m)
            # Tuple-marginal product measures have tuple variates:
            x isa Tuple ? (@test length(x) == length(d)) : (@test size(x) == size(d))
            @test insupport(m, x) != false
        end
    end

    @testset "Dirac" begin
        d = Distributions.Dirac(4.2)
        m = @inferred asmeasure(d)
        @test m === MeasureBase.Dirac(4.2)
        @test_throws ArgumentError AsMeasure{typeof(d)}(d)
        @test @inferred(Distributions.Distribution(m)) === d
    end

    @testset "products" begin
        test_conversion(product_distribution(Weibull.([0.7, 1.1, 1.3])), ProductMeasure)
        test_conversion(product_distribution(Poisson.([0.7, 1.4])), ProductMeasure)

        if isdefined(Distributions, :ProductDistribution)
            test_conversion(product_distribution(Weibull(0.7), Exponential(1.3)), ProductMeasure)
        end
    end

    @testset "reshaped" begin
        test_conversion(reshape(MvNormal([0.7, 0.9], [1.4 0.5; 0.5 1.1]), 1, 2), PushforwardMeasure)
        test_conversion(reshape(product_distribution(Weibull.([0.7, 1.1, 1.3, 0.9, 1.2, 0.8])), 2, 3), PushforwardMeasure)
    end

    @testset "mixtures" begin
        test_conversion(MixtureModel([Normal(-1.0, 1.0), Normal(2.0, 3.0)], [0.3, 0.7]), SuperpositionMeasure)
        test_conversion(MixtureModel([Normal(-2.0, 1.0), Normal(0.0, 2.0), Normal(3.0, 1.0)], [0.2, 0.5, 0.3]), SuperpositionMeasure)
        test_conversion(MixtureModel([Exponential(0.3), Weibull(2.0, 1.0)], [0.4, 0.6]), SuperpositionMeasure)
        test_conversion(MixtureModel([MvNormal([0.0, 0.0], I(2)), MvNormal([2.0, 2.0], 2 * I(2))], [0.3, 0.7]), SuperpositionMeasure)
        test_conversion(UnivariateGMM([-1.0, 2.0], [1.0, 0.5], Categorical([0.4, 0.6])), SuperpositionMeasure)

        d = MixtureModel([Normal(-1.0, 1.0), Normal(2.0, 3.0)], [0.3, 0.7])
        m = asmeasure(d)
        @test probs(convert(Distribution, m)) ≈ probs(d)
        @test massof(m) ≈ 1
        @test massof(asmeasure(Normal())) == 1
        @test mean(rand(stblrng(), Float64, m^1000)) ≈ mean(d) atol = 0.3

        # Hand-built superpositions of weighted probability measures behave
        # like mixtures:
        m2 = 0.3 * asmeasure(Normal(-1.0, 1.0)) + 0.7 * asmeasure(Normal(2.0, 3.0))
        for x in (rand(stblrng(), d) for _ in 1:10)
            @test logdensityof(m2, x) ≈ logpdf(d, x)
        end
    end

    @testset "standard distributions" begin
        @test StandardUniform === StandardDist{Uniform}
        @test StandardNormal === StandardDist{Normal}
        @test StandardUniform{0} === StandardDist{Uniform,0}
        @test StandardNormal{1} === StandardDist{Normal,1}

        for (D, B) in [
            (Uniform, StdUniform()),
            (Exponential, StdExponential()),
            (Logistic, StdLogistic()),
            (Normal, StdNormal()),
        ]
            @test @inferred(asmeasure(StandardDist{D}())) === B
            @test @inferred(asmeasure(StandardDist{D}(3))) == B^3
            @test @inferred(asmeasure(StandardDist{D}(2, 3))) == B^(2, 3)

            @test @inferred(Distributions.Distribution(B)) === StandardDist{D}()
            @test @inferred(convert(Distribution, B)) === StandardDist{D}()
            @test @inferred(Distributions.Distribution(B^3)) == StandardDist{D}(3)
            @test @inferred(convert(Distribution, B^(2, 3))) == StandardDist{D}(2, 3)

            d = StandardDist{D}(3)
            x = rand(stblrng(), d)
            @test logdensityof(asmeasure(d), x) ≈ logpdf(d, x)
        end
    end
end
