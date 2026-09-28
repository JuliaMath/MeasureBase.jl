# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random
using StableRNGs: StableRNG
using OneTwoMany: firstarg, secondarg

using MeasureBase
using MeasureBase: StdExponential, StdLogistic, StdNormal, StdUniform
using MeasureBase: mcombine, mbind, productmeasure, transport_to, pushfwd
using StructArrays: StructArray
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
        Xab = vcat(randn(2, 3), rand(2, 3))
        @test logdensities(mpab, Xab) ≈ [logdensityof(mpab, x) for x in eachcol(Xab)]

        # Products of array-variate marginals combine into flat streams:
        N2 = StdNormal()^2
        pv = mcombine(vcat, productmeasure([N2, N2]), productmeasure([N2]))
        x6 = randn(6)
        @test logdensityof(pv, x6) ≈ logdensityof(StdNormal()^6, x6)
        @test length(rand(pv)) == 6
        X6 = randn(6, 3)
        @test logdensities(pv, X6) ≈ [logdensityof(StdNormal()^6, x) for x in eachcol(X6)]
        u6 = transport_to(StdUniform()^6, pv)(x6)
        @test u6 ≈ transport_to(StdUniform()^6, StdNormal()^6)(x6)
        @test transport_to(pv, StdUniform()^6)(u6) ≈ x6

        # Merging products with overlapping names is an error:
        pm = mcombine(merge, productmeasure((a = StdNormal(),)), productmeasure((b = StdUniform(),)))
        @test logdensityof(pm, (a = 0.3, b = 0.5)) ≈ logdensityof(StdNormal(), 0.3)
        @test_throws ArgumentError mcombine(merge, productmeasure((a = StdNormal(), b = StdNormal())), productmeasure((b = StdUniform(),)))
    end

    @testset "batched densities of merge combinations" begin
        N, E = StdNormal(), StdExponential()
        Bm = mbind(p -> productmeasure((c = pushfwd(Mul(abs(p.a) + 1), N), d = E^length(p.b))), productmeasure((a = N, b = E^2)), merge)
        Cm = mcombine(merge, Bm, productmeasure((e = N,)))
        @test Cm isa MeasureBase.CombinedMeasure
        for M in (Bm, Cm)
            xs = [rand(stblrng(), M) for _ in 1:4]
            @test logdensities(M, StructArray(xs)) ≈ logdensityof.(Ref(M), xs)
        end
    end

    @testset "nested vcat combinations" begin
        T = productmeasure((StdNormal(), StdExponential()))
        C = mcombine(vcat, mcombine(vcat, mcombine(vcat, T, StdNormal()), StdNormal())^2, StdNormal())
        x = rand(stblrng(), C)
        @test length(x) == 9
        ℓ_ref = sum(logdensityof(StdNormal(), x[i]) for i in (1, 3, 4, 5, 7, 8, 9)) + sum(logdensityof(StdExponential(), x[i]) for i in (2, 6))
        @test logdensityof(C, x) ≈ ℓ_ref
        @test logdensities(C, hcat(x, x)) ≈ [ℓ_ref, ℓ_ref]
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
