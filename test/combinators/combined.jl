# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random
using StableRNGs: StableRNG
using OneTwoMany: firstarg, secondarg

using MeasureBase
using MeasureBase: StdExponential, StdLogistic, StdNormal, StdUniform
using MeasureBase: mcombine, productmeasure, transport_to, pushfwd
using AffineMaps: Mul

@testset "mcombine" begin
    stblrng() = StableRNG(789990641)

    α = StdExponential()
    β = StdLogistic()

    @testset "combination shortcuts" begin
        @test mcombine(firstarg, α, β) === α
        @test mcombine(secondarg, α, β) === β
        @test mcombine(tuple, α, β) == productmeasure((α, β))
        @test mcombine(vcat, StdNormal()^2, StdNormal()^1) == StdNormal()^3
        @test mcombine(vcat, productmeasure([α, α]), productmeasure([α])) ==
              productmeasure([α, α, α])
        @test mcombine(
            merge,
            productmeasure((a = α,)),
            productmeasure((b = β,)),
        ) == productmeasure((a = α, b = β))
        @test mcombine(tuple, MeasureBase.Dirac(1), MeasureBase.Dirac(2)) ==
              MeasureBase.Dirac((1, 2))
    end

    @testset "mcombine of products" begin
        p1 = productmeasure([pushfwd(Mul(1.0), StdNormal()), pushfwd(Mul(2.0), StdNormal())])
        p2 = productmeasure([pushfwd(Mul(3.0), StdNormal())])
        p12 = mcombine(vcat, p1, p2)
        @test p12 isa MeasureBase.ProductMeasure && length(MeasureBase.marginals(p12)) == 3
        x3 = randn(3)
        @test logdensityof(p12, x3) ≈ logdensityof(p1, x3[1:2]) + logdensityof(p2, x3[3:3])

        # Products of different or abstract marginal types stay combined measures:
        mab = mcombine(vcat, StdNormal()^2, StdUniform()^2)
        @test mab isa MeasureBase.CombinedMeasure
        pa = productmeasure(AbstractMeasure[StdNormal(), StdNormal()])
        pb = productmeasure(AbstractMeasure[StdUniform(), StdUniform()])
        mpab = mcombine(vcat, pa, pb)
        @test mpab isa MeasureBase.CombinedMeasure
        x = vcat(randn(2), rand(2))
        @test logdensityof(mpab, x) ≈ logdensityof(mab, x)
        @test_throws ArgumentError logdensities(mpab, vcat(randn(2, 3), rand(2, 3)))
    end

    @testset "CombinedMeasure" begin
        μ = mcombine(Pair, α, β)
        @test μ isa MeasureBase.CombinedMeasure

        ab = rand(stblrng(), Float64, μ)
        @test ab isa Pair
        @test logdensityof(μ, ab) ≈ logdensityof(α, ab.first) + logdensityof(β, ab.second)

        @test MeasureBase.getdof(μ) == 2
        @test MeasureBase.fast_dof(μ) == 2
        @test MeasureBase.insupport(μ, ab) isa MeasureBase.NoFastInsupport

        y = transport_to(StdUniform()^2, μ)(ab)
        @test y isa AbstractVector{<:Real} && length(y) == 2
        ab_reco = transport_to(μ, StdUniform()^2)(y)
        @test ab_reco.first ≈ ab.first && ab_reco.second ≈ ab.second
    end
end
