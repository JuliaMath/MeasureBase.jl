using Test

using MeasureBase
using MeasureBase: StdNormal, Half, SpikeMixture, logdensities
using MeasureBase: pushfwd, superpose, PushfwdRootMeasure
import MeasureBase
using Distributions: Normal, Poisson, Binomial, Dirichlet, MvNormal, logpdf, pdf
using AffineMaps: Mul, Add
using InverseFunctions: InverseFunctions
using IntervalSets
using LinearAlgebra: I
using StaticArrays: SVector
using Static: static

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

    @testset "powers of primitive roots" begin
        x = [0.1, -0.2, 0.3]
        @test logdensity_rel(StdNormal()^static(3), StdNormal()^3, x) == 0
        @test logdensity_rel(StdNormal()^3, StdNormal()^static(3), SVector(x...)) == 0
        @test logdensity_rel(asmeasure(MvNormal(zeros(3), I(3))), StdNormal()^static(3), SVector(x...)) ≈ 0 atol = 1e-14
        @test logdensity_rel(P^2, StdNormal()^2, [1.0, 2.0]) == Inf
        @test logdensity_rel(StdNormal()^2, P^2, [1.0, 2.0]) == -Inf
        @test logdensity_rel(StdNormal()^2, P^2, [1.5, 2.0]) == Inf
        p2 = pushfwd(Mul(2.0), StdNormal(), PushfwdRootMeasure())
        @test logdensity_rel(p2^3, StdNormal()^3, x) ≈ sum(logpdf.(Normal(0, 2), x) .- logpdf.(Normal(), x))
    end

    @testset "superpositions of continuous and atomic measures" begin
        s = superpose(P^2, StdNormal()^2)
        @test logdensityof(s, [1.0, 2.0]) ≈ logpdf(Poisson(2.0), 1) + logpdf(Poisson(2.0), 2)
        @test logdensityof(s, [1.5, 2.0]) ≈ logpdf(Normal(), 1.5) + logpdf(Normal(), 2.0)

        v = superpose(StdNormal()^2, Dirac([0.0, 0.0]))
        ℓ_cont = logpdf(Normal(), 0.5) + logpdf(Normal(), 0.0)
        @test logdensityof(v, [0.0, 0.0]) == 0
        @test logdensityof(v, [0.5, 0.0]) ≈ ℓ_cont
        @test logdensities(v, [0.0 0.5; 0.0 0.0]) ≈ [0, ℓ_cont]
        @test logdensityof(superpose(StdNormal(), Dirac(0.0)), 0.0) == 0

        d3 = superpose(StdNormal()^3, Dirac(0.5)^3)
        @test logdensityof(d3, fill(0.5, 3)) == 0
        @test logdensityof(d3, [0.5, 0.5, 0.4]) ≈ sum(logpdf.(Normal(), [0.5, 0.5, 0.4]))

        # Both components have the same density:
        sv = superpose(StdNormal()^static(2), asmeasure(MvNormal(zeros(2), I(2))))
        @test logdensity_rel(sv, StdNormal()^static(2), SVector(0.1, 0.2)) ≈ log(2)
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
        @test logdensity_rel(pp^3, StdNormal()^3, [2.0, 4.0, 6.0]) == Inf
        @test logdensity_rel(pushfwd(Mul(2.0), Dirac(0.0) + Dirac(1.0)), Dirac(0.0) + Dirac(2.0), 2.0) == 0

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

    @testset "pushforwards of mixed measures" begin
        ps = pushfwd(Mul(2.0), StdNormal() + Dirac(0.5))
        @test rootmeasure(ps) == MeasureBase.LebesgueBase() + Dirac(1.0)
        @test logdensityof(ps, 1.0) == 0
        @test logdensityof(ps, 0.5) ≈ logpdf(Normal(0, 2), 0.5)
        @test logdensity_rel(ps, pushfwd(Mul(2.0), StdNormal()) + Dirac(1.0), 1.0) == 0
        @test logdensity_rel(ps, pushfwd(Mul(2.0), StdNormal()) + Dirac(1.0), 0.3) ≈ 0 atol = 1e-14

        n = MeasureBase.normalize(0.25 * StdNormal() + 0.75 * Dirac(0.5))
        @test logdensityof(pushfwd(Add(1.0), n), 1.5) ≈ log(0.75)
        @test logdensityof(pushfwd(Add(1.0), n), 1.0) ≈ log(0.25) + logpdf(Normal(), 0.0)

        # Mixed roots are pushed forward as a whole:
        sm = SpikeMixture(StdNormal(), 0.7)
        @test logdensityof(pushfwd(Mul(2.0), sm), 0.0) ≈ log(0.3)
        @test logdensityof(pushfwd(Add(1.0), sm), 1.0) ≈ log(0.3)
        @test logdensityof(pushfwd(Add(1.0), sm), 0.0) ≈ log(0.7) + logpdf(Normal(), -1.0)
    end
end
