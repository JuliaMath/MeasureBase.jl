# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random
using StableRNGs: StableRNG
using AffineMaps: Mul

using MeasureBase
using MeasureBase: StdExponential, StdNormal, StdUniform
using MeasureBase: mbind, mkernel, bindkernel, boundmeasure
using MeasureBase: pushfwd, productmeasure, transport_to, transportmeasure, localmeasure

@testset "bind" begin
    stblrng() = StableRNG(789990641)

    f_β(σ) = pushfwd(Mul(σ + 0.5), StdNormal())
    α = StdExponential()

    @testset "monadic bind" begin
        μ = mbind(f_β, α)
        @test μ isa MeasureBase.Bind
        @test boundmeasure(μ) === α
        @test bindkernel(μ) isa MeasureBase.MKernel
        @test mbind(bindkernel(μ), α) == μ
        @test mbind(f_β)(α) == μ

        a = rand(stblrng(), Float64, α)
        b = rand(copy(stblrng()), Float64, μ)  # not comparable directly, just smoke:
        @test rand(stblrng(), Float64, μ) isa Real

        @test MeasureBase.insupport(μ, 0.4) isa MeasureBase.NoFastInsupport
        @test MeasureBase.getdof(μ) isa MeasureBase.NoDOF
        @test_throws ArgumentError basemeasure(μ)
        @test_throws ArgumentError MeasureBase.rootmeasure(μ)
    end

    @testset "mbind with tuple and Pair" begin
        for f_c in (tuple, Pair)
            μ = mbind(f_β, α, f_c)
            ab = rand(stblrng(), Float64, μ)
            a, b = f_c === tuple ? ab : (ab.first, ab.second)
            @test logdensityof(μ, ab) ≈ logdensityof(α, a) + logdensityof(f_β(a), b)

            tpm = transportmeasure(μ, ab)
            @test logdensityof(tpm, ab) ≈ logdensityof(μ, ab)
            @test localmeasure(μ, ab) == tpm
        end
    end

    @testset "mbind with vcat" begin
        αv = StdExponential()^1
        f_βv(a) = pushfwd(Mul(a[1] + 0.5), StdNormal())^2
        μ = mbind(f_βv, αv, vcat)

        xy = rand(stblrng(), Float64, μ)
        @test xy isa AbstractVector{<:Real} && length(xy) == 3
        a, b = xy[1:1], xy[2:3]
        @test logdensityof(μ, xy) ≈ logdensityof(αv, a) + logdensityof(f_βv(a), b)

        # Transport to and from a standard measure, dof of μ is not fast-computable:
        y = transport_to(StdUniform()^3, μ)(xy)
        @test y isa AbstractVector{<:Real} && length(y) == 3
        @test all(u -> 0 <= u <= 1, y)
        xy_reco = transport_to(μ, StdUniform()^3)(y)
        @test xy_reco ≈ xy
    end

    @testset "mbind with merge" begin
        αnt = productmeasure((position = StdNormal(),))
        f_βnt(a) = productmeasure((
            noise = pushfwd(Mul(abs(a.position) + 0.5), StdExponential()),
        ))
        μ = mbind(f_βnt, αnt, merge)

        x = rand(stblrng(), Float64, μ)
        @test x isa NamedTuple{(:position, :noise)}
        @test logdensityof(μ, x) ≈
              logdensityof(αnt, (position = x.position,)) +
              logdensityof(f_βnt(x), (noise = x.noise,))

        y = transport_to(StdUniform()^2, μ)(x)
        @test y isa AbstractVector{<:Real} && length(y) == 2
        x_reco = transport_to(μ, StdUniform()^2)(y)
        @test x_reco.position ≈ x.position && x_reco.noise ≈ x.noise
    end

    @testset "mbind with Dirac" begin
        @test mbind(f_β, MeasureBase.Dirac(1.5)) == asmeasure(f_β(1.5))
        μ = mbind(f_β, MeasureBase.Dirac(1.5), tuple)
        ab = rand(stblrng(), Float64, μ)
        @test ab[1] == 1.5
    end
end
