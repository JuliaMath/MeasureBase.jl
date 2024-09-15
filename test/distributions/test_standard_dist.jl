# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using DistributionMeasures
using Test

using Random, Statistics, LinearAlgebra
using Distributions, PDMats
using StableRNGs
import ForwardDiff, ChainRulesTestUtils


@testset "standard_dist" begin
    stblrng() = StableRNG(789990641)

    for (D, sz, dref) in [
        (Uniform, (), Uniform()),
        (Uniform, (5,), product_distribution(fill(Uniform(0.0, 1.0), 5))),
        (Uniform, (2, 3), reshape(product_distribution(fill(Uniform(0.0, 1.0), 6)), 2, 3)),
        (Normal, (), Normal()),
        (Normal, (), Normal(0., 1.0)),
        (Normal, (5,), MvNormal(Diagonal(fill(1.0, 5)))),
        (Normal, (2, 3), reshape(MvNormal(Diagonal(fill(1.0, 6))), 2, 3)),
        (Exponential, (), Exponential()),
        (Exponential, (5,), product_distribution(fill(Exponential(1.0), 5))),
        (Exponential, (2, 3), reshape(product_distribution(fill(Exponential(1.0), 6)), 2, 3)),
    ]
        @testset "StandardDist{$D}($(join(sz,",")))" begin
            N = length(sz)

            @test @inferred(StandardDist{D}(sz...)) isa StandardDist{D}
            @test @inferred(StandardDist{D}(sz...)) isa StandardDist{D}
            @test @inferred(size(StandardDist{D}(sz...))) == size(dref)
            @test @inferred(size(StandardDist{D}(sz...))) == size(dref)

            d = StandardDist{D}(sz...)

            if size(d) == ()
                @test @inferred(DistributionMeasures.nonstddist(d)) == dref
            end

            @test @inferred(length(d)) == length(dref)
            @test @inferred(size(d)) == size(dref)

            @test @inferred(eltype(typeof(d))) == eltype(typeof(dref))
            @test @inferred(eltype(d)) == eltype(dref)

            @test @inferred(Distributions.params(d)) == ()
            @test @inferred(partype(d)) == partype(dref)

            for f in [minimum, maximum, mean, median, mode, modes, var, std, skewness, kurtosis, location, scale, entropy]
                supported_by_dref = try f(dref); true catch MethodError; false; end
                if supported_by_dref
                    @test @inferred(f(d)) ≈ f(dref)
                end
            end

            for x in [rand(dref) for i in 1:10]
                ref_gradlogpdf = try
                    gradlogpdf(dref, x)
                catch MethodError
                    ForwardDiff.gradient(x -> logpdf(dref, x), x)
                end
                @test @inferred(gradlogpdf(d, x)) ≈ ref_gradlogpdf
                @test @inferred(logpdf(d, x)) ≈ logpdf(dref, x)
                @test @inferred(pdf(d, x)) ≈ pdf(dref, x)
            end

            if size(d) == ()
                for x in [minimum(dref), quantile(dref, 1//3), quantile(dref, 1//2), quantile(dref, 2//3), maximum(dref)]
                    for f in [logpdf, pdf, gradlogpdf, logcdf, cdf, logccdf, ccdf]
                        @test @inferred(f(d, x)) ≈ f(dref, x)
                    end
                end

                for x in [0, 1//3, 1//2, 2//3, 1]
                    for f in [quantile, cquantile]
                        @test @inferred(f(d, x)) ≈ f(dref, x)
                    end
                end

                for x in log.([0, 1//3, 1//2, 2//3, 1])
                    for f in [invlogcdf, invlogccdf]
                        @test @inferred(f(d, x)) ≈ f(dref, x)
                    end
                end

                for p in [0.0, 0.25, 0.75, 1.0]
                    @test @inferred(quantile(d, p)) == quantile(dref, p)
                    @test @inferred(cquantile(d, p)) == cquantile(dref, p)
                end

                for t in [-3, 0, 3]
                    @test isapprox(@inferred(mgf(d, t)), mgf(dref, t), rtol = 1e-5)
                    @test isapprox(@inferred(cf(d, t)), cf(dref, t), rtol = 1e-5)
                end

                @test @inferred(truncated(d, quantile(dref, 1//3), quantile(dref, 2//3))) == truncated(dref, quantile(dref, 1//3), quantile(dref, 2//3))

                @test @inferred(product_distribution(fill(d, 3))) == StandardDist{typeof(d)}(3)
                @test @inferred(product_distribution(fill(d, 3, 4))) == StandardDist{typeof(d)}(3, 4)
            end

            if length(size(d)) == 1
                @test @inferred(convert(Distributions.Product, d)) isa Distributions.Product
                d_as_prod = convert(Distributions.Product, d)
                @test d_as_prod.v == fill(StandardDist{D}(), size(d)...)
            end

            @test @inferred(rand(stblrng(), d)) == rand(stblrng(), d)
            @test @inferred(rand(stblrng(), d, 5)) == rand(stblrng(), d, 5)

            @test @inferred(rand(stblrng(), d)) == rand(stblrng(), dref)
            @test @inferred(rand(stblrng(), d, 5)) == rand(stblrng(), dref, 5)
            @test @inferred(rand!(stblrng(), d, zeros(size(d)...))) == rand!(stblrng(), dref, zeros(size(dref)...))
            if length(size(d)) == 1
                @test @inferred(rand!(stblrng(), d, zeros(size(d)..., 5))) == rand!(stblrng(), dref, zeros(size(dref)..., 5))
            end
        end
    end

    @testset "StandardDist{Normal}()" begin
        # TODO: Add @inferred
        d = StandardDist{Normal}(4)
        d_uv = StandardDist{Normal}()
        dref = MvNormal(Diagonal(fill(1.0, 4)))
        @test (MvNormal(d)) == dref
        @test (Base.convert(MvNormal, d)) == dref
    end
end
