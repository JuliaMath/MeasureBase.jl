using Test
using MeasureBase
using MeasureBase: insupport as measure_insupport

using DensityInterface: logdensityof
import Distributions: insupport as dist_insupport
using Distributions: Normal, Exponential, Logistic, Uniform

@testset "standard" begin
    for (m, d) in [
        (StdUniform(), Uniform()),
        (StdExponential(), Exponential()),
        (StdLogistic(), Logistic()),
        (StdNormal(), Normal()),
    ]
        @testset "$(nameof(typeof(m)))" begin
            for x in [-Inf, -1.2, -1, 0, 0.0, 1 // 2, 0.5, 1, 1.0, 2, 2.3]
                @test @inferred(logdensityof(m, x)) ≈ logdensityof(d, x)
                @test @inferred(logdensityof(m, x)) ≈ logdensity_rel(m, rootmeasure(m), x)
                @test measure_insupport(m, x) ≈ dist_insupport(d, x)
                if measure_insupport(m, x)
                    @test @inferred(MeasureBase.unsafe_logdensityof(m, x)) ≈
                          logdensityof(d, x)
                end
            end
        end
    end
end

using LogExpFunctions: log1pexp
using IrrationalConstants: log2π

@testset "Standard measure logdensityof" begin
    # StdNormal
    @test logdensityof(StdNormal(), 0.0) ≈ -log2π / 2
    @test logdensityof(StdNormal(), 1.0) ≈ (-1 - log2π) / 2

    # StdUniform
    @test logdensityof(StdUniform(), 0.5) == 0.0
    @test logdensityof(StdUniform(), 1.5) == -Inf
    @test logdensityof(StdUniform(), -0.5) == -Inf

    # StdLogistic
    x = 1.0
    @test logdensityof(StdLogistic(), x) ≈ -abs(x) - 2 * log1pexp(-abs(x))

    # StdExponential
    @test logdensityof(StdExponential(), 1.0) ≈ -1.0
    @test logdensityof(StdExponential(), -1.0) == -Inf
end
