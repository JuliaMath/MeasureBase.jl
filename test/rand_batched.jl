# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random, Statistics
using StableRNGs: StableRNG
using StructArrays: StructArray
using ArraysOfArrays: flatview

using MeasureBase
using MeasureBase: GenContext, get_precision
using MeasureBase.HeterogeneousComputing: AbstractComputeUnit
using MeasureBase: StdUniform, StdExponential, StdLogistic, StdNormal, Dirac
using MeasureBase: weightedmeasure, superpose, mcombine, mbind, productmeasure, pushfwd, SpikeMixture
using MeasureBase: rand_impl, batched_rand_impl, asmeasure, logdensities
using Distributions: Normal, MvNormal, Weibull, Poisson, logpdf
using JLArrays: JLArrays, JLArray, JLBackend

struct NoBatchRandMeasure <: AbstractMeasure end
MeasureBase.rand_impl(ctx::GenContext, ::NoBatchRandMeasure) = 3 * rand(MeasureBase.get_rng(ctx), MeasureBase.get_precision(ctx))

struct NoRandMeasure <: AbstractMeasure end

@testset "batched rand" begin
    stblrng() = StableRNG(789)
    ctx() = GenContext{Float64}(stblrng())

    @testset "single variates as batches with zero batch dimensions" begin
        @test @inferred(batched_rand_impl(ctx(), StdNormal(), ())) isa Float64
        @test batched_rand_impl(ctx(), StdNormal(), ()) == rand(stblrng(), StdNormal())
        @test @inferred(batched_rand_impl(ctx(), StdExponential(), ())) isa Float64
        @test @inferred(batched_rand_impl(ctx(), StdLogistic(), ())) isa Float64
        @test @inferred(batched_rand_impl(ctx(), StdUniform(), ())) isa Float64
        @test @inferred(batched_rand_impl(ctx(), StdNormal()^3, ())) isa Vector{Float64}
        @test batched_rand_impl(ctx(), StdNormal()^3, ()) == rand(stblrng(), StdNormal()^3)
        @test batched_rand_impl(ctx(), (StdNormal()^2)^3, ()) == flatview(rand(stblrng(), (StdNormal()^2)^3))
        @test batched_rand_impl(ctx(), Dirac(2.5), ()) === 2.5
        @test batched_rand_impl(ctx(), Dirac([1.0, 2.0]), ()) == [1.0, 2.0]
        @test batched_rand_impl(ctx(), weightedmeasure(0.3, StdNormal()), ()) isa Float64
        @test batched_rand_impl(ctx(), pushfwd(exp, StdNormal()), ()) isa Float64
        @test batched_rand_impl(ctx(), pushfwd(Base.BroadcastFunction(exp), StdNormal()^2), ()) isa Vector{Float64}
        @test batched_rand_impl(ctx(), SpikeMixture(StdNormal(), 0.5), ()) isa Float64
        @test batched_rand_impl(ctx(), superpose(StdNormal(), StdUniform()), ()) isa Float64
        @test rand(ctx(), NoBatchRandMeasure()) isa Float64
        @test batched_rand_impl(ctx(), NoBatchRandMeasure(), (4,)) isa Vector{Float64}
        @test_throws ArgumentError rand(NoRandMeasure())
        @test_throws ArgumentError batched_rand_impl(ctx(), NoRandMeasure(), ())
        @test_throws ArgumentError batched_rand_impl(ctx(), NoRandMeasure(), (3,))
    end

    @testset "integer draws" begin
        # Integer variates stay integers, also as sizes in binds:
        P = asmeasure(Poisson(3.0))
        @test rand(stblrng(), Float32, P) isa Int
        @test eltype(rand(stblrng(), Float32, P^3)) == Int
        m = mbind(n -> StdNormal()^n, P, tuple)
        n, x = rand(stblrng(), m)
        @test n isa Int && length(x) == n
        @test logdensityof(m, (n, x)) ≈ logpdf(Poisson(3.0), n) + sum(logdensityof.(StdNormal(), x))
        @test length(rand(stblrng(), mbind(n -> StdNormal()^n, P))) isa Int
    end

    @testset "structured batches" begin
        Pt = productmeasure((StdNormal(), StdExponential()^2))
        Xt = batched_rand_impl(ctx(), Pt, (5,))
        @test Xt isa Tuple && size(Xt[1]) == (5,) && size(Xt[2]) == (2, 5)
        xt = batched_rand_impl(ctx(), Pt, ())
        @test xt isa Tuple{Float64,Vector{Float64}} && xt == rand(stblrng(), Pt)
        XPt = rand(stblrng(), Pt^5)
        @test XPt isa StructArray && length(XPt) == 5
        @test XPt[2] isa Tuple{Float64,<:AbstractVector{Float64}} && length(XPt[2][2]) == 2
        @test logdensityof(Pt^5, collect(XPt)) ≈ sum(logdensityof.(Ref(Pt), XPt))
        Pn = productmeasure((a = StdNormal(), b = StdExponential()^2))
        Xn = batched_rand_impl(ctx(), Pn, (2, 3))
        @test Xn isa NamedTuple{(:a, :b)} && size(Xn.a) == (2, 3) && size(Xn.b) == (2, 2, 3)
        XPn = rand(stblrng(), Pn^4)
        @test XPn isa StructArray && XPn[1] isa NamedTuple{(:a, :b)}

        mm = mcombine(merge, productmeasure((a = StdNormal(),)), productmeasure((b = StdUniform()^2,)))
        Xm = batched_rand_impl(ctx(), mm, (3,))
        @test Xm isa NamedTuple{(:a, :b)} && size(Xm.a) == (3,) && size(Xm.b) == (2, 3)
        @test batched_rand_impl(ctx(), mm, ()) == rand(stblrng(), mm)
        mt = mcombine(tuple, StdNormal(), StdUniform()^2)
        Xtt = batched_rand_impl(ctx(), mt, (3,))
        @test Xtt isa Tuple && size(Xtt[1]) == (3,) && size(Xtt[2]) == (2, 3)
        @test batched_rand_impl(ctx(), mt, ()) == rand(stblrng(), mt)
    end

    @testset "powers and value-dependent sizes" begin
        f_β(a) = StdNormal()^length(a)
        μb = mbind(f_β, StdUniform()^1, vcat)
        x = rand(stblrng(), μb^3)
        @test x isa Vector{Vector{Float64}} && length(x) == 3 && all(length.(x) .== 2)
        X = batched_rand_impl(ctx(), μb^3, (4,))
        @test size(X) == (2, 3, 4)
        @test batched_rand_impl(ctx(), μb, ()) == rand(stblrng(), μb)
        @test size(batched_rand_impl(ctx(), (StdNormal()^2)^3, (4, 5))) == (2, 3, 4, 5)
        @test size(batched_rand_impl(ctx(), StdNormal()^(2, 3), (4,))) == (2, 3, 4)
    end

    @testset "moments of batches" begin
        n = 20_000
        X = batched_rand_impl(ctx(), superpose(StdNormal(), Dirac(3.0)), (n,))
        @test isapprox(mean(X), 1.5, atol = 0.05)
        Xs = batched_rand_impl(ctx(), SpikeMixture(StdNormal()^2, 0.5), (n,))
        @test size(Xs) == (2, n) && isapprox(mean(Xs .== 0), 0.5, atol = 0.02)
        Xp = batched_rand_impl(ctx(), pushfwd(Base.BroadcastFunction(exp), StdNormal()^2), (n,))
        @test isapprox(mean(log.(Xp)), 0.0, atol = 0.03)
        Xd = batched_rand_impl(ctx(), MeasureBase.AsMeasure{Normal{Float64}}(Normal(1.0, 2.0)), (n,))
        @test Xd isa Vector{Float64} && isapprox(mean(Xd), 1.0, atol = 0.05)
        xd = batched_rand_impl(ctx(), MeasureBase.AsMeasure{Normal{Float64}}(Normal(1.0, 2.0)), ())
        @test xd isa Float64
        Xmv = batched_rand_impl(ctx(), MeasureBase.AsMeasure{typeof(MvNormal([1.0, 2.0], [1.0, 0.5]))}(MvNormal([1.0, 2.0], [1.0, 0.5])), (n,))
        @test size(Xmv) == (2, n) && isapprox(vec(mean(Xmv, dims = 2)), [1.0, 2.0], atol = 0.05)
    end

    @testset "variates on device compute units" begin
        JLArrays.allowscalar(false)
        dctx = GenContext{Float32}(AbstractComputeUnit(JLBackend()))
        Random.seed!(MeasureBase.get_rng(dctx), 5678)
        n = 10_000
        for (μ, m) in (
            (StdNormal(), 0.0), (StdUniform(), 0.5), (StdExponential(), 1.0), (StdLogistic(), 0.0),
            (StdNormal()^3, 0.0), ((StdNormal()^2)^3, 0.0), (weightedmeasure(0.3, StdExponential()^2), 1.0),
            (mcombine(vcat, StdNormal()^2, StdUniform()^1), 1 / 6), (pushfwd(Base.BroadcastFunction(exp), StdNormal()^2), exp(0.5)),
            (superpose(StdNormal(), StdUniform()), 0.25), (SpikeMixture(StdExponential(), 0.5), 0.5),
            (asmeasure(Normal(0.3, 1.7)), 0.3), (asmeasure(Weibull(1.4, 0.9)), mean(Weibull(1.4, 0.9))),
        )
            X = batched_rand_impl(dctx, μ, (n,))
            @test X isa JLArray{Float32}
            @test all(isfinite, Array(logdensities(μ, X)))
            @test isapprox(mean(Array(X)), m, atol = 0.05)
            @test get_precision(rand(dctx, μ)) == Float32
        end
        Pt = productmeasure((StdNormal(), StdExponential()^2))
        Xt = batched_rand_impl(dctx, Pt, (n,))
        @test Xt isa Tuple{JLArray{Float32,1},JLArray{Float32,2}} && size(Xt[2]) == (2, n)
        @test all(isfinite, Array(logdensities(Pt, Xt)))
        @test isapprox(mean(Array(Xt[2])), 1.0, atol = 0.05)
    end
end
