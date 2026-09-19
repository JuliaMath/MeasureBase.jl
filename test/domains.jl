# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: combinesets, setcartprod, setcartpower
using MeasureBase: CombinedSet, CartesianProduct, CartesianPower, ImplicitDomain
using MeasureBase: maybe_in, mdomain, mcombine, mbind, pushfwd
using MeasureBase: StdNormal, StdUniform, StdExponential
using MeasureBase: ℝ, ℤ
using OneTwoMany: firstarg, secondarg
using AffineMaps: Mul

@testset "domains" begin
    @testset "combinesets" begin
        @test combinesets(firstarg, ℝ, ℤ) === ℝ
        @test combinesets(secondarg, ℝ, ℤ) === ℤ

        s_tuple = combinesets(tuple, ℝ, ℤ)
        @test s_tuple isa CartesianProduct
        @test (1.5, 2) ∈ s_tuple
        @test !((1.5, 2.5) ∈ s_tuple)

        pv = setcartprod([ℝ, ℝ])
        sv = combinesets(vcat, pv, pv)
        @test sv isa CartesianProduct
        @test [1.0, 2.0, 3.0, 4.0] ∈ sv

        snt = combinesets(merge, setcartprod((a = ℝ,)), setcartprod((b = ℤ,)))
        @test snt isa CartesianProduct
        @test (a = 1.5, b = 2) ∈ snt

        # One-dimensional powers of equal singleton base sets concatenate:
        spw = combinesets(vcat, setcartpower(ℝ, (2,)), setcartpower(ℝ, (3,)))
        @test spw isa CartesianPower
        @test [1.0, 2.0, 3.0, 4.0, 5.0] ∈ spw
        @test combinesets(vcat, setcartpower(ℝ, (2,)), setcartpower(ℤ, (3,))) isa
              CombinedSet

        # No specific representation available:
        sc = combinesets(vcat, ℝ, setcartpower(ℝ, (2,)))
        @test sc isa CombinedSet
        @test maybe_in([1.0, 2.0, 3.0], sc)
        @test !isempty(sc)
        @test_throws ArgumentError [1.0, 2.0, 3.0] ∈ sc

        # Implicit domains combine into the implicit domain of the
        # combined measure:
        sid = combinesets(
            vcat,
            ImplicitDomain(StdNormal()^2),
            ImplicitDomain(StdUniform()^1),
        )
        @test sid isa ImplicitDomain
        @test maybe_in([1.0, 2.0, 3.0], sid)
    end

    @testset "mdomain of combined measures" begin
        f_β(a) = pushfwd(Mul(a[1] + 0.5), StdNormal())^2
        μc = mcombine(vcat, StdNormal()^2, mbind(f_β, StdExponential()^1, vcat))
        @test μc isa MeasureBase.CombinedMeasure
        @test mdomain(μc) isa ImplicitDomain
        @test maybe_in(rand(Float64, μc), mdomain(μc))
    end
end
