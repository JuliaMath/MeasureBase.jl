using Test

using MeasureBase
using MeasureBase: StdNormal, Half, logdensities
using MeasureBase: pushfwd, superpose, productmeasure, PushfwdRootMeasure
import MeasureBase
using Distributions: Normal, Poisson, Binomial, Beta, Dirichlet, MvNormal, logpdf, pdf
using AffineMaps: Mul, Add
using InverseFunctions: InverseFunctions, setinverse
using IntervalSets
using LinearAlgebra: I
using StaticArrays: SVector
using Static: static
using SpecialFunctions: besseli, besselix

struct _SimplexWarp{T}
    r::T
end
(f::_SimplexWarp)(x) = [x[1]^f.r, 1 - x[1]^f.r]
InverseFunctions.inverse(f::_SimplexWarp) = _SimplexWarp(inv(f.r))

# The von Mises distribution on the unit circle in the plane, with a
# density relative to arc length, a reference measure without a chart:
struct _VonMisesCircle <: MeasureBase.AbstractMeasure
    κ::Float64
end
MeasureBase.basemeasure(::_VonMisesCircle) = MeasureBase.OpaqueRoot(:unit_circle, (2,), static(1))
MeasureBase.logdensity_def(m::_VonMisesCircle, x) = m.κ * (x[1] - 1) - log(2π * besselix(0, m.κ))
MeasureBase.insupportof(::_VonMisesCircle, x) = x[1]^2 + x[2]^2 ≈ 1
MeasureBase.mspace_ndims(::Type{<:_VonMisesCircle}) = 1

struct _OnSphere <: MeasureBase.AbstractMeasure end
MeasureBase.basemeasure(::_OnSphere) = MeasureBase.OpaqueRoot(:unit_sphere, (3,), static(2))
MeasureBase.logdensity_def(::_OnSphere, x) = -log(4π)
MeasureBase.insupportof(::_OnSphere, x) = sum(abs2, x) ≈ 1

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

    @testset "products of primitive roots" begin
        s = superpose(productmeasure((StdNormal(), StdNormal())), Dirac((0.0, 1.0)))
        @test logdensityof(s, (0.0, 1.0)) == 0
        @test logdensityof(s, (0.5, 1.0)) ≈ logpdf(Normal(), 0.5) + logpdf(Normal(), 1.0)
        c = productmeasure([Counting(1:3), Counting(0:4)])
        @test logdensity_rel(c, StdNormal()^2, [1.0, 2.0]) == Inf
        @test logdensity_rel(c, P^2, [1.0, 2.0]) ≈ -(logpdf(Poisson(2.0), 1) + logpdf(Poisson(2.0), 2))
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

        spike = 2 * Dirac(0.0)
        slab = asmeasure(Beta(0.5, 1))
        mixture = spike + slab
        @test densityof(mixture, 0.0) ≈ 2
        @test density_rel(mixture, spike, 0.0) ≈ 1
        @test logdensities(mixture, [0.0, 0.25]) ≈ [log(2), 0]

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

    @testset "opaque roots" begin
        m1, m2 = _VonMisesCircle(1.0), _VonMisesCircle(2.5)
        x = [cos(0.3), sin(0.3)]
        vm(κ, θ) = exp(κ * cos(θ)) / (2π * besseli(0, κ))
        @test logdensityof(m1, x) ≈ log(vm(1.0, 0.3))
        @test logdensity_rel(m1, m2, x) ≈ log(vm(1.0, 0.3) / vm(2.5, 0.3))
        @test logdensity_rel(m1^2, m2^2, [x x]) ≈ 2 * log(vm(1.0, 0.3) / vm(2.5, 0.3))
        @test logdensityof(m1 + m2, x) ≈ log(vm(1.0, 0.3) + vm(2.5, 0.3))
        @test getdof(m1) == 1 && MeasureBase.fast_dof(m1) == 1 && getdof(m1^3) == 3
        @test getdof(_OnSphere()) == 2
        @test getdof(MeasureBase.OpaqueRoot(:unknown, (3, 3))) isa MeasureBase.NoDOF

        # Other roots don't relate:
        @test_throws ArgumentError logdensity_rel(m1, StdNormal()^2, x)
        @test_throws ArgumentError logdensity_rel(StdNormal()^2, m1, x)
        @test_throws ArgumentError logdensity_rel(m1, MeasureBase.Lebesgue()^2, x)
        @test_throws ArgumentError logdensity_rel(m1, MeasureBase.OpaqueRoot(:ellipse, (2,), static(1)), x)
        @test_throws ArgumentError logdensityof(m1 + StdNormal()^2, x)

        # Maps of the ambient space have no volume element for opaque roots,
        # the root can only be pushed forward:
        f = setinverse(x -> 2x, y -> y / 2)
        @test_throws ArgumentError pushfwd(f, m1)
        @test_throws ArgumentError pushfwd(f, m1^3)
        @test_throws ArgumentError pushfwd(f, productmeasure((m1, StdNormal())))
        p = pushfwd(f, m1, PushfwdRootMeasure())
        @test logdensityof(p, 2x) ≈ logdensityof(m1, x)
        @test_throws ArgumentError logdensity_rel(p, m1, x)
    end

    @testset "pushforwards of mixed measures" begin
        to_uniform = transport_to(StdUniform(), StdNormal())
        mixed = StdNormal() + 2 * Dirac(0.0)
        uniform_spike = pushfwd(to_uniform, mixed)
        @test density_rel(uniform_spike, Dirac(0.5), 0.5) ≈ 2
        @test densityof(uniform_spike, 0.75) ≈ 1
        @test density_rel(pushfwd(to_uniform, MeasureBase.normalize(mixed)), Dirac(0.5), 0.5) ≈ 2 / 3
        @test density_rel(pushfwd(to_uniform, P), Dirac(0.5), 0.5) ≈ pdf(Poisson(2.0), 0)

        ps = pushfwd(Mul(2.0), StdNormal() + Dirac(0.5))
        @test rootmeasure(ps) == MeasureBase.MeasureJoin((MeasureBase.LebesgueBase(), Dirac(1.0)))
        @test logdensityof(ps, 1.0) == 0
        @test logdensityof(ps, 0.5) ≈ logpdf(Normal(0, 2), 0.5)
        @test logdensity_rel(ps, pushfwd(Mul(2.0), StdNormal()) + Dirac(1.0), 1.0) == 0
        @test logdensity_rel(ps, pushfwd(Mul(2.0), StdNormal()) + Dirac(1.0), 0.3) ≈ 0 atol = 1e-14

        n = MeasureBase.normalize(0.25 * StdNormal() + 0.75 * Dirac(0.5))
        @test logdensityof(pushfwd(Add(1.0), n), 1.5) ≈ log(0.75)
        @test logdensityof(pushfwd(Add(1.0), n), 1.0) ≈ log(0.25) + logpdf(Normal(), 0.0)

        sm = 0.7 * StdNormal() + 0.3 * Dirac(0.0)
        @test logdensityof(pushfwd(Mul(2.0), sm), 0.0) ≈ log(0.3)
        @test logdensityof(pushfwd(Add(1.0), sm), 1.0) ≈ log(0.3)
        @test logdensityof(pushfwd(Add(1.0), sm), 0.0) ≈ log(0.7) + logpdf(Normal(), -1.0)
        f = setinverse(x -> (2 * x[1], x[2] + 1), y -> (y[1] / 2, y[2] - 1))
        mp = pushfwd(f, productmeasure((StdNormal(), P)))
        @test logdensityof(mp, (1.0, 3.0)) ≈ logpdf(Normal(), 0.5) + logpdf(Poisson(2.0), 2)

        line = mcombine(vcat, Dirac(1.0), StdNormal())
        for scale in ([2.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 2.0])
            image = pushfwd(Mul(scale), line)
            y = scale * [1.0, 0.5]
            @test logdensityof(rootmeasure(image), y) == 0
            @test logdensityof(image, y) ≈ logpdf(Normal(), 0.5)
        end
    end
end
