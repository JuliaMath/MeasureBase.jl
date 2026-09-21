# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random, Statistics
using StableRNGs: StableRNG
using Static: static
using ArraysOfArrays: flatview
using AffineMaps: Add

using MeasureBase
using MeasureBase: GenContext
using MeasureBase: StdUniform, StdExponential, StdLogistic, StdNormal, Dirac, Lebesgue
using MeasureBase: weightedmeasure, superpose, mcombine, mbind, productmeasure, pushfwd, testvalue
using MeasureBase: rand_impl, batched_rand_impl, massof, isnormalized

@testset "rand" begin
    stblrng() = StableRNG(789990641)

    @testset "generative contexts" begin
        @test rand(stblrng(), StdNormal()) == rand(stblrng(), StdNormal())
        @test rand(stblrng(), StdNormal()) == rand(GenContext{Float64}(stblrng()), StdNormal())
        @test rand(stblrng(), StdNormal()^3) == rand(stblrng(), StdNormal()^3)
        @test @inferred(rand(stblrng(), Float32, StdNormal())) isa Float32
        @test @inferred(rand(Float32, StdNormal()^3)) isa Vector{Float32}
        @test @inferred(rand(GenContext{Float32}(stblrng()), StdUniform()^(2, 3))) isa Matrix{Float32}
        @test @inferred(rand(StdExponential())) isa Float64
        @test_throws ArgumentError rand(Lebesgue())
    end

    @testset "layout of power variates" begin
        x = rand(StdNormal()^(2, 3))
        @test x isa Matrix{Float64} && size(x) == (2, 3)
        xs = rand(StdNormal()^static(3))
        @test xs isa AbstractVector{Float64} && length(xs) == 3
        xn = rand((StdNormal()^2)^3)
        @test xn isa AbstractVector && length(xn) == 3 && all(x -> length(x) == 2, xn)
        @test size(flatview(xn)) == (2, 3)
        @test logdensityof((StdNormal()^2)^3, xn) ≈ logdensityof(StdNormal()^6, vec(flatview(xn)))
        @test batched_rand_impl(GenContext{Float64}(stblrng()), StdNormal()^2, (4, 5)) isa Array{Float64,3}
        @test size(batched_rand_impl(GenContext{Float64}(stblrng()), (StdNormal()^2)^3, (4,))) == (2, 3, 4)
    end

    @testset "test values" begin
        @test testvalue(mcombine(vcat, StdNormal()^2, StdUniform()^3)) == [0.0, 0.0, 0.5, 0.5, 0.5]
        @test rand(MeasureBase.ConstantRNG(), Float64, StdUniform()^3) == fill(0.5, 3)
        @test rand(MeasureBase.ConstantRNG(), Float32, StdLogistic()^2) == zeros(Float32, 2)
        @test testvalue(StdNormal()) == 0
        @test testvalue(Float32, StdUniform()) === 0.5f0
        @test testvalue(StdExponential()^3) == ones(3)
        @test testvalue((StdLogistic()^2)^2) == [zeros(2), zeros(2)]
        @test testvalue(productmeasure((a = StdNormal(), b = StdUniform()^2))) == (a = 0.0, b = [0.5, 0.5])
    end

    @testset "distribution of variates" begin
        n = 10^5
        for (μ, m, v) in [
            (StdUniform(), 0.5, 1 / 12),
            (StdExponential(), 1.0, 1.0),
            (StdLogistic(), 0.0, π^2 / 3),
            (StdNormal(), 0.0, 1.0),
            (weightedmeasure(0.3, StdNormal()), 0.0, 1.0),
            (pushfwd(exp, StdNormal()), exp(0.5), (exp(1) - 1) * exp(1)),
            (MeasureBase.Half(StdNormal()), sqrt(2 / π), 1 - 2 / π),
            (SpikeMixture(Dirac(1.0), 0.25), 0.25, 0.25 * 0.75),
        ]
            X = rand(stblrng(), μ^n)
            @test isapprox(mean(X), m, atol = 5 * sqrt(v / n) + 1e-3)
            @test isapprox(var(X), v, rtol = 0.05)
            xs = [rand_impl(GenContext{Float64}(stblrng()), μ) for _ in 1:20]
            @test all(x -> insupport(μ, x) != false, xs)
        end

        mix = superpose(weightedmeasure(log(0.3), Dirac(0.0)), weightedmeasure(log(0.7), Dirac(1.0)))
        @test isapprox(mean(rand(stblrng(), mix^n)), 0.7, atol = 0.01)
        @test isapprox(mean([rand(mix) for _ in 1:n]), 0.7, atol = 0.01)
        mixn = superpose(weightedmeasure(log(0.5), StdNormal()), weightedmeasure(log(0.5), pushfwd(Add(4.0), StdNormal())))
        @test isapprox(mean(rand(stblrng(), mixn^n)), 2.0, atol = 0.02)
    end

    @testset "structural measures" begin
        P = MeasureBase.ProductMeasure([weightedmeasure(log(i), StdNormal()) for i in 1:3])
        @test @inferred(rand(stblrng(), P)) isa Vector{Float64}
        XP = rand(stblrng(), P^100)
        @test size(flatview(XP)) == (3, 100)
        @test logdensities(P, XP) ≈ [logdensityof(P, x) for x in XP]

        Pt = productmeasure((StdNormal(), StdUniform()^2))
        xt = rand(stblrng(), Pt)
        @test xt isa Tuple{Float64,Vector{Float64}}
        Xt = rand(stblrng(), Pt^5)
        @test Xt isa AbstractVector && length(Xt) == 5

        mc = mcombine(vcat, StdNormal()^2, StdUniform()^3)
        xc = rand(stblrng(), mc)
        @test xc isa Vector{Float64} && length(xc) == 5 && all(0 .<= xc[3:5] .<= 1)
        Xc = rand(stblrng(), mc^50)
        @test size(flatview(Xc)) == (5, 50)
        @test logdensities(mc, Xc) ≈ [logdensityof(mc, x) for x in Xc]

        f_β(a) = StdNormal()^length(a)
        μb = mbind(f_β, StdUniform()^2, vcat)
        xb = rand(stblrng(), μb)
        @test xb isa AbstractVector && length(xb) == 4
        Xb = rand(stblrng(), μb^3)
        @test Xb isa AbstractVector && length(Xb) == 3 && all(x -> length(x) == 4, Xb)

        ctx = GenContext{Float64}(stblrng())
        spd = superpose(weightedmeasure(log(0.5), Dirac([1.0, 2.0])), weightedmeasure(log(0.5), Dirac([3.0, 4.0])))
        Xspd = batched_rand_impl(ctx, spd, (6,))
        @test size(Xspd) == (2, 6) && all(c -> c == [1.0, 2.0] || c == [3.0, 4.0], eachcol(Xspd))
        Xsm = batched_rand_impl(ctx, SpikeMixture(StdNormal()^3, 0.5), (4,))
        @test size(Xsm) == (3, 4) && all(c -> all(iszero, c) || !any(iszero, c), eachcol(Xsm))
        @test massof(StdNormal()) == 1 && massof(StdUniform()^3) == 1 && massof(weightedmeasure(log(2.0), StdNormal()^2)) ≈ 2
        @test isnormalized(StdNormal()) && isnormalized(StdExponential()^(2, 2)) && !isnormalized(2.0 * StdNormal())
        @test !isnormalized(Lebesgue())

        d = Dirac([1.0, 2.0])
        @test rand(d^2) == [[1.0, 2.0], [1.0, 2.0]]
        @test batched_rand_impl(GenContext{Float64}(stblrng()), d, (2,)) == [1.0 1.0; 2.0 2.0]
    end
end
