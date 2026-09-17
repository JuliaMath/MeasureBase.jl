# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: logdensities, logdensity_def, StdNormal, StdUniform, Dirac, Lebesgue, LebesgueBase, superpose, weightedmeasure
using ArraysOfArrays: VectorOfSimilarVectors, sliced, flatview
using StaticArrays: SVector, @SVector, @SMatrix
using Static: static
using IrrationalConstants: log2π
import JLArrays
using JLArrays: JLArray

stdnormal_ld(x) = -(x^2 + log2π) / 2

@testset "logdensities" begin
    @testset "scalar variates" begin
        X = randn(10)
        @test @inferred(logdensities(StdNormal(), X)) ≈ stdnormal_ld.(X)
        Xm = randn(2, 3)
        @test logdensities(StdNormal(), Xm) ≈ stdnormal_ld.(Xm)
    end

    @testset "powers with nested variates" begin
        m3 = StdNormal()^3
        X = [randn(3) for _ in 1:10]
        @test @inferred(logdensities(m3, X)) ≈ [sum(stdnormal_ld, x) for x in X]
        @test only(logdensities(m3, [X[1]])) ≈ logdensityof(m3, X[1])

        m23 = StdNormal()^(2, 3)
        X23 = [randn(2, 3) for _ in 1:5]
        @test logdensities(m23, X23) ≈ [sum(stdnormal_ld, x) for x in X23]

        mpp = (StdNormal()^(2, 3))^4
        Xpp = [[randn(2, 3) for _ in 1:4] for _ in 1:6]
        @test logdensities(mpp, Xpp) ≈ [sum(x -> sum(stdnormal_ld, x), xs) for xs in Xpp]
    end

    @testset "powers with flat variate storage" begin
        m3 = StdNormal()^3
        X = VectorOfSimilarVectors(randn(3, 10))
        @test @inferred(logdensities(m3, X)) ≈
              vec(sum(stdnormal_ld.(flatview(X)), dims = 1))

        # Power structure may be stored flattened out within each point:
        mpp = (StdNormal()^(2, 3))^4
        Xpp = sliced(randn(2, 3, 4, 7), 3)
        @test logdensities(mpp, Xpp) ≈ [sum(stdnormal_ld, x) for x in Xpp]
    end

    @testset "non-scalar-variate fallback" begin
        mprod = productmeasure((StdUniform(), StdNormal()))
        X = [(rand(), randn()) for _ in 1:5]
        @test logdensities(mprod, X) ≈ logdensityof.(Ref(mprod), X)
    end

    @testset "unknown variate size" begin
        mix = superpose(weightedmeasure(log(0.3), StdNormal()), weightedmeasure(log(0.7), StdUniform()))
        X = randn(4, 5)
        @test logdensities(mix, X) ≈ logdensityof.(Ref(mix), X)
        @test logdensityof(mix^4, X[:, 1]) ≈ sum(logdensityof.(Ref(mix), X[:, 1]))
        @test logdensities(mix^4, sliced(X, 1)) ≈ vec(sum(logdensityof.(Ref(mix), X), dims = 1))
    end

    @testset "size mismatch" begin
        @test_throws ArgumentError logdensities(StdNormal()^3, [randn(3), randn(2)])
        @test_throws ArgumentError logdensities(
            StdNormal()^3,
            VectorOfSimilarVectors(randn(2, 5)),
        )
    end

    @testset "flat batch storage and array-variate bases" begin
        m3 = StdNormal()^3
        Xf = randn(3, 10)
        @test @inferred(logdensities(m3, Xf)) ≈ vec(sum(stdnormal_ld.(Xf), dims = 1))
        x = randn(3)
        @test @inferred(logdensities(m3, x)) ≈ sum(stdnormal_ld, x)

        mpp = (StdNormal()^(2, 3))^4
        Xpp_flat = randn(2, 3, 4, 7)
        @test @inferred(logdensities(mpp, Xpp_flat)) ≈ vec(sum(stdnormal_ld.(Xpp_flat), dims = (1, 2, 3)))
        Xpp_nested = sliced(sliced(Xpp_flat, 2), 1)
        @test logdensities(mpp, Xpp_nested) ≈ logdensities(mpp, Xpp_flat)
        xpp = randn(2, 3, 4)
        @test @inferred(logdensityof(mpp, xpp)) ≈ logdensityof(mpp, [xpp[:, :, i] for i in 1:4])
        @test logdensityof(mpp, sliced(xpp, 2)) ≈ logdensityof(mpp, xpp)

        mvec = Dirac([1.0, 2.0])^3
        @test @inferred(logdensities(mvec, [fill([1.0, 2.0], 3) for _ in 1:2])) == [0.0, 0.0]
    end

    @testset "static variates" begin
        m3 = StdNormal()^static(3)
        xs = @SVector randn(3)
        f(x) = logdensityof(m3, x)
        @test @inferred(f(xs)) ≈ sum(stdnormal_ld, xs)
        @test @allocated(f(xs)) == 0
        g(x) = logdensityof(StdNormal()^3, x)
        xd = randn(3)
        @test @inferred(g(xd)) ≈ sum(stdnormal_ld, xd)
        @test @allocated(g(xd)) == 0
        Xs = @SMatrix randn(3, 4)
        @test @inferred(logdensities(m3, Xs)) ≈ vec(sum(stdnormal_ld.(Xs), dims = 1))
        @test logdensities(m3, Xs) isa SVector{4}
        @test @inferred(logdensityof(StdNormal()^static(0), SVector{0,Float64}())) == 0
        @test @inferred(logdensityof(StdNormal()^0, Float64[])) == 0
    end

    @testset "powers of primitive measures" begin
        @test @inferred(logdensity_def(Lebesgue()^3, randn(3))) == 0
        @test @inferred(logdensity_def(LebesgueBase()^(2, 2), randn(2, 2))) == 0
        @test @inferred(logdensityof(Lebesgue()^3, randn(3))) == 0
    end

    @testset "GPU array semantics" begin
        JLArrays.allowscalar(false)

        X = JLArray(randn(100))
        ld = logdensities(StdNormal(), X)
        @test ld isa JLArray
        @test Array(ld) ≈ stdnormal_ld.(Array(X))

        Xb = VectorOfSimilarVectors(JLArray(randn(3, 50)))
        ldb = logdensities(StdNormal()^3, Xb)
        @test ldb isa JLArray
        @test Array(ldb) ≈ vec(sum(stdnormal_ld.(Array(flatview(Xb))), dims = 1))

        xj = JLArray(randn(10))
        @test logdensityof(StdNormal()^10, xj) ≈ logdensityof(StdNormal()^10, Array(xj))
    end
end
