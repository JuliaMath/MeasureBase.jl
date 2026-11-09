using Test

using MeasureBase
using MeasureBase: StdNormal, Half
using MeasureBase: pushfwd, superpose, PushfwdRootMeasure
import MeasureBase
using Distributions: Normal, Poisson, Binomial, Dirichlet, logpdf, pdf
using AffineMaps: Mul, Add
using InverseFunctions: InverseFunctions
using IntervalSets

struct _SimplexWarp{T}
    r::T
end
(f::_SimplexWarp)(x) = [x[1]^f.r, 1 - x[1]^f.r]
InverseFunctions.inverse(f::_SimplexWarp) = _SimplexWarp(inv(f.r))

@testset "root_relations.jl" begin
    P = asmeasure(Poisson(2.0))

    @testset "pushed-forward Lebesgue roots" begin
        p2 = pushfwd(Mul(2.0), StdNormal(), PushfwdRootMeasure())
        p3 = pushfwd(Mul(3.0), StdNormal(), PushfwdRootMeasure())
        for y in (-1.5, 0.3, 1.0)
            @test logdensity_rel(p2, p3, y) ≈ logpdf(Normal(0, 2), y) - logpdf(Normal(0, 3), y)
            @test logdensity_rel(p2, StdNormal(), y) ≈ logpdf(Normal(0, 2), y) - logpdf(Normal(), y)
        end
    end

    @testset "pushed-forward counting roots" begin
        pp = pushfwd(Mul(2.0), P)
        @test logdensityof(pp, 2.0) ≈ logpdf(Poisson(2.0), 1)
        @test logdensityof(pp, 3.0) == -Inf
        @test insupportof(pp, 2.0) && !insupportof(pp, 3.0)
        @test logdensityof(superpose(pp, StdNormal()), 2.5) ≈ logpdf(Normal(), 2.5)
        @test logdensityof(superpose(pp, StdNormal()), 2.0) ≈ logpdf(Poisson(2.0), 1)
        @test logdensity_rel(pp, Dirac(2.0), 2.0) ≈ logpdf(Poisson(2.0), 1)
        @test logdensity_rel(pp, StdNormal(), 2.0) == Inf
        @test logdensity_rel(StdNormal(), pp, 2.0) == -Inf

        P3 = asmeasure(Poisson(3.0))
        q2, q3 = pushfwd(Mul(2.0), P), pushfwd(Mul(3.0), P3)
        @test logdensity_rel(q2, q3, 3.0) == -Inf
        @test logdensity_rel(q2, q3, 6.0) ≈ logpdf(Poisson(2.0), 3) - logpdf(Poisson(3.0), 2)
        s = superpose(q2, q3)
        @test logdensityof(s, 3.0) ≈ logpdf(Poisson(3.0), 1)
        @test logdensity_rel(s, CountingBase(), 6.0) ≈ log(pdf(Poisson(2.0), 3) + pdf(Poisson(3.0), 2))
    end

    @testset "distinct pushed-forward references don't cancel" begin
        d = asmeasure(Dirichlet([1.0, 1.0]))
        m, n = pushfwd(_SimplexWarp(2.0), d), pushfwd(_SimplexWarp(3.0), d)
        y = [0.25, 0.75]
        @test_throws ArgumentError logdensity_rel(m, n, y)
        @test logdensity_rel(m, pushfwd(_SimplexWarp(2.0), d), y) == 0
    end

    @testset "pushed-forward atoms" begin
        b = pushfwd(Add(-1.0), asmeasure(Binomial(2, 0.5)))
        @test_throws ArgumentError Half(b)
        @test_throws ArgumentError truncate(b, 0.0 .. 1.0)
    end
end
