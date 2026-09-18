# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: logdensities, logdensity_def, StdNormal, StdUniform, StdExponential, StdLogistic, Dirac, Lebesgue, LebesgueBase, superpose, weightedmeasure, mcombine, productmeasure
using ArraysOfArrays: VectorOfSimilarVectors, sliced, flatview
using StaticArrays: SVector, @SVector, @SMatrix
using Static: static
using IrrationalConstants: log2π
import JLArrays
using JLArrays: JLArray

stdnormal_ld(x) = -(x^2 + log2π) / 2

# A measure with array variates of a known flat size at the type level:
struct VecTestMeasure{T} <: AbstractMeasure
    s::T
end
MeasureBase.mspace_elsize(::VecTestMeasure) = (2,)
MeasureBase.mspace_flatsize(::VecTestMeasure) = (2,)
MeasureBase.mspace_flatsize(::Type{<:VecTestMeasure}) = (2,)
MeasureBase.basemeasure(::VecTestMeasure) = LebesgueBase()^2
MeasureBase.insupport(::VecTestMeasure, x) = true
MeasureBase.logdensityof_impl(m::VecTestMeasure, x) = -sum(abs2, x) / (2 * m.s)

# `@allocated` at top level reports a boxed result on Julia 1.10:
_allocated(f::F, args::Vararg{Any,N}) where {F,N} = @allocated f(args...)

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
        @test @inferred(logdensityof(mix^0, Float64[])) == 0
        @test logdensities(mix^0, [Float64[], Float64[]]) == [0.0, 0.0]
    end

    @testset "flat and nested variate forms agree" begin
        for (μ, x_flat) in (
            ((StdNormal()^3)^2, randn(3, 2)),
            ((StdNormal()^(2, 3))^4, randn(2, 3, 4)),
        )
            x_nested = sliced(x_flat, length(MeasureBase.mspace_flatsize(MeasureBase.pwr_base(μ))))
            @test logdensityof(μ, x_flat) ≈ logdensityof(μ, x_nested)
            @test MeasureBase.checked_arg(μ, x_flat) === x_flat
            @test MeasureBase.checked_arg(μ, x_nested) === x_nested
            @test_throws ArgumentError MeasureBase.checked_arg(μ, randn(7))
        end
    end

    @testset "size mismatch" begin
        @test_throws ArgumentError logdensityof(StdNormal()^3, randn(3, 4))
        @test_throws ArgumentError logdensityof(StdNormal()^(2, 3), randn(2, 3, 1))
        @test_throws ArgumentError logdensityof(StdNormal()^3, randn(4))
        @test_throws ArgumentError logdensityof(StdNormal()^3, 1.0)
        @test_throws ArgumentError logdensities(StdNormal(), VectorOfSimilarVectors(randn(3, 5)))
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
        @test @inferred(logdensityof(m3, xs)) ≈ sum(stdnormal_ld, xs)
        @test _allocated(logdensityof, m3, xs) == 0
        xd = randn(3)
        @test @inferred(logdensityof(StdNormal()^3, xd)) ≈ sum(stdnormal_ld, xd)
        @test _allocated(logdensityof, StdNormal()^3, xd) == 0
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

    @testset "structural batched kernels" begin
        w = weightedmeasure(log(0.3), StdNormal()^3)
        X = randn(3, 10)
        @test @inferred(logdensities(w, X)) ≈ [logdensityof(w, x) for x in eachcol(X)]
        xw = randn(3)
        @test @inferred(logdensityof(w, xw)) ≈ log(0.3) + sum(stdnormal_ld, xw)
        @test _allocated(logdensityof, w, xw) == 0

        ms = [weightedmeasure(log(i), StdNormal()) for i in 1:4]
        prod4 = productmeasure(ms)
        @test @inferred(MeasureBase.mspace_flatsize(prod4)) == (4,)
        @test @inferred(MeasureBase.mspace_elsize(prod4)) == (4,)
        xp = randn(4)
        @test @inferred(logdensityof(prod4, xp)) ≈ sum(log(i) + stdnormal_ld(xp[i]) for i in 1:4)
        @test _allocated(logdensityof, prod4, xp) == 0
        Xp = randn(4, 7)
        @test @inferred(logdensities(prod4, Xp)) ≈ [logdensityof(prod4, x) for x in eachcol(Xp)]
        @test logdensities(prod4, sliced(Xp, 1)) ≈ logdensities(prod4, Xp)
        @test @inferred(logdensityof(prod4^2, randn(4, 2))) isa Float64
        Xpp = randn(4, 2, 5)
        @test logdensities(prod4^2, Xpp) ≈ [logdensityof(prod4^2, Xpp[:, :, i]) for i in 1:5]

        ms2 = reshape([weightedmeasure(log(i), StdUniform()) for i in 1:6], 2, 3)
        prod23 = productmeasure(ms2)
        @test @inferred(MeasureBase.mspace_flatsize(prod23)) == (2, 3)
        x23 = rand(2, 3)
        @test @inferred(logdensityof(prod23, x23)) ≈ sum(log(i) for i in 1:6)
        @test logdensities(prod23, rand(2, 3, 4)) ≈ fill(sum(log(i) for i in 1:6), 4)

        mvec = productmeasure([StdNormal()^2, StdNormal()^2])
        @test @inferred(MeasureBase.mspace_flatsize(mvec)) isa MeasureBase.NoMSpaceElementSize
    end

    @testset "batched with-rest for combined measures" begin
        m = mcombine(vcat, StdNormal()^2, StdUniform()^3)
        @test @inferred(MeasureBase.mspace_flatsize(m)) == (5,)
        x = vcat(randn(2), rand(3))
        @test @inferred(logdensityof(m, x)) ≈ sum(stdnormal_ld, x[1:2])
        X = vcat(randn(2, 6), rand(3, 6))
        @test @inferred(logdensities(m, X)) ≈ [logdensityof(m, x) for x in eachcol(X)]
        @test logdensities(m, sliced(X, 1)) ≈ logdensities(m, X)
        ℓ, A_rest = MeasureBase.batched_logdensityof_with_rest(StdNormal()^2, X, ())
        @test ℓ ≈ vec(sum(stdnormal_ld.(X[1:2, :]), dims = 1))
        @test size(A_rest) == (3, 6)
        ℓ2, A_rest2 = MeasureBase.batched_logdensityof_with_rest(StdNormal(), X, (2,))
        @test MeasureBase._materialize(ℓ2) ≈ stdnormal_ld.(X[1:2, :]) && size(A_rest2) == (3, 6)
        @test_throws ArgumentError logdensities(m, vcat(X, rand(1, 6)))

        m3 = mcombine(vcat, StdNormal(), mcombine(vcat, StdExponential()^2, StdLogistic()))
        @test @inferred(MeasureBase.mspace_flatsize(m3)) == (4,)
        X3 = vcat(randn(1, 5), rand(2, 5), randn(1, 5))
        @test @inferred(logdensities(m3, X3)) ≈ [logdensityof(m3, x) for x in eachcol(X3)]
    end

    @testset "static sizes in combined measures" begin
        m = mcombine(vcat, StdNormal()^static(2), StdUniform()^3)
        @test @inferred(MeasureBase.mspace_flatsize(m)) == (5,)
        x = vcat(randn(2), rand(3))
        @test @inferred(logdensityof(m, x)) ≈ sum(stdnormal_ld, x[1:2])
        X = vcat(randn(2, 4), rand(3, 4))
        @test @inferred(logdensities(m, X)) ≈ [logdensityof(m, x) for x in eachcol(X)]
        ms = mcombine(vcat, StdNormal()^static(2), StdUniform()^static(3))
        @test @inferred(MeasureBase.mspace_flatsize(ms)) == MeasureBase.mspace_flatsize(StdNormal()^static(5))
        @test @inferred(logdensityof(ms, SVector{5}(x))) ≈ logdensityof(m, x)
        @test @inferred(logdensityof(ms, x)) ≈ logdensityof(m, x)
        @test logdensities(ms, X) ≈ logdensities(m, X)
        @test_throws ArgumentError MeasureBase.batched_logdensityof_with_rest(StdNormal(), zeros(0, 4), ())
    end

    @testset "array products of array-variate marginals" begin
        p = productmeasure([VecTestMeasure(1.0), VecTestMeasure(2.0), VecTestMeasure(0.5)])
        @test @inferred(MeasureBase.mspace_flatsize(p)) == (2, 3)
        xs = [randn(2) for _ in 1:3]
        X = stack(xs)
        ℓ = sum(map(logdensityof, MeasureBase.marginals(p), xs))
        @test @inferred(logdensityof(p, xs)) ≈ ℓ
        @test @inferred(logdensityof(p, X)) ≈ ℓ
        A = randn(2, 3, 5)
        @test @inferred(logdensities(p, A)) ≈ [logdensityof(p, A[:, :, i]) for i in 1:5]
        @test logdensities(p, sliced(A, Val(2))) ≈ logdensities(p, A)
        @test logdensities(p, zeros(2, 3, 0)) == Float64[]
        @test_throws ArgumentError logdensityof(p, randn(2, 2))
        @test_throws ArgumentError logdensities(p, randn(2, 2, 5))

        # Powers with static axes have a flat size at the type level:
        pp = MeasureBase.ProductMeasure([StdNormal()^static(2), StdNormal()^static(2)])
        @test @inferred(MeasureBase.mspace_flatsize(pp)) == (2, 2)
        Xp = randn(2, 2)
        @test @inferred(logdensityof(pp, Xp)) ≈ sum(stdnormal_ld, Xp)
        A3 = randn(2, 2, 3)
        @test @inferred(logdensities(pp, A3)) ≈ [sum(stdnormal_ld, A3[:, :, i]) for i in 1:3]
    end

    @testset "products with mixed marginal types" begin
        pa = productmeasure(AbstractMeasure[StdNormal(), StdUniform()])
        X = vcat(randn(1, 4), rand(1, 4))
        xs = [X[:, i] for i in 1:4]
        @test logdensities(pa, xs) ≈ [logdensityof(pa, x) for x in xs]
        @test logdensities(pa, sliced(X, Val(1))) ≈ [logdensityof(pa, x) for x in xs]
        @test_throws ArgumentError logdensityof(pa, 0.5)
        @test_throws ArgumentError logdensityof(productmeasure((StdNormal(), StdUniform())), 0.5)
        @test_throws ArgumentError logdensityof(pa, X[:, 1:1])
    end

    @testset "out-of-support and empty batches of structural kernels" begin
        w = weightedmeasure(log(0.3), StdUniform()^2)
        X = [0.5 -0.5 0.5; 0.5 0.5 1.5]
        @test logdensities(w, X) == [log(0.3), -Inf, -Inf]
        @test logdensities(w, zeros(2, 0)) == Float64[]
        pu = productmeasure([weightedmeasure(log(i), StdUniform()) for i in 1:2])
        @test logdensities(pu, X) == [log(2), -Inf, -Inf]
        mc = mcombine(vcat, StdUniform()^1, StdExponential()^1)
        @test logdensities(mc, [0.5 -0.5 0.5; 0.5 0.5 -1.0]) == [-0.5, -Inf, -Inf]
        @test logdensities(mc, zeros(2, 0)) == Float64[]
        @test MeasureBase.logdensity_def(pu, [0.5, 0.5]) ≈ log(2)
    end

    @testset "GPU array semantics for structural kernels" begin
        JLArrays.allowscalar(false)
        ms = JLArray([weightedmeasure(log(i), StdNormal()) for i in 1:4])
        prodj = productmeasure(ms)
        Xj = JLArray(randn(4, 7))
        ldj = logdensities(prodj, Xj)
        @test ldj isa JLArray
        @test Array(ldj) ≈ logdensities(productmeasure(Array(ms)), Array(Xj))
        mc = mcombine(vcat, StdNormal()^2, StdUniform()^3)
        Xc = JLArray(vcat(randn(2, 6), rand(3, 6)))
        ldc = logdensities(mc, Xc)
        @test ldc isa JLArray
        @test Array(ldc) ≈ logdensities(mc, Array(Xc))
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
