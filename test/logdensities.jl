# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: logdensities, StdNormal, StdUniform
using ArraysOfArrays: VectorOfSimilarVectors, nestedview, flatview
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
        Xpp = nestedview(randn(2, 3, 4, 7), 3)
        @test logdensities(mpp, Xpp) ≈ [sum(stdnormal_ld, x) for x in Xpp]
    end

    @testset "non-scalar-variate fallback" begin
        mprod = productmeasure((StdUniform(), StdNormal()))
        X = [(rand(), randn()) for _ in 1:5]
        @test logdensities(mprod, X) ≈ logdensityof.(Ref(mprod), X)
    end

    @testset "size mismatch" begin
        @test_throws ArgumentError logdensities(StdNormal()^3, [randn(3), randn(2)])
        @test_throws ArgumentError logdensities(
            StdNormal()^3,
            VectorOfSimilarVectors(randn(2, 5)),
        )
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
