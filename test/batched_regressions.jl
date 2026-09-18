# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Regression tests for the batched-first review findings.

using Test

using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, Dirac, Lebesgue, GenContext
using MeasureBase: productmeasure, pushfwd, mcombine, mbind, weightedmeasure, insupport, testvalue
using MeasureBase: batched_rand_impl, batched_transport_to_std_with_rest, batched_transport_from_std_with_rest
using MeasureBase.InverseFunctions: inverse
using ArraysOfArrays: VectorOfVectors, nestedview, flatview, sliced
using StaticArrays: SVector
using Static: static
using Distributions: MvNormal, logpdf
using LinearAlgebra: I
using JLArrays

struct UnknownRankMeasure <: AbstractMeasure end
MeasureBase.basemeasure(::UnknownRankMeasure) = Lebesgue()
MeasureBase.logdensity_def(::UnknownRankMeasure, x) = -sum(abs, x)
MeasureBase.insupport(::UnknownRankMeasure, x) = true

# A parameterized function object (not unwrapped into struct array columns):
struct Scale <: Function
    s::Float64
end
(f::Scale)(x) = f.s * x
MeasureBase.InverseFunctions.inverse(f::Scale) = Scale(inv(f.s))
MeasureBase.ChangesOfVariables.with_logabsdet_jacobian(f::Scale, x) = (f(x), log(abs(f.s)))

@testset "batched regressions" begin
    @testset "products of function-wrapper marginals" begin
        pm = productmeasure([pushfwd(Scale(s), StdExponential()) for s in 0.1:0.2:0.9])
        x = rand(pm)
        @test logdensityof(pm, x) ≈ sum(logdensityof.(MeasureBase.marginals(pm), x))
        X = rand(pm^4)
        @test logdensities(pm, X) ≈ [logdensityof(pm, X[j]) for j in 1:4]
        f = transport_to(StdUniform()^5, pm)
        @test flatview(inverse(f).(f.(sliced(flatview(X), Val(1))))) ≈ flatview(X)
        JLArrays.allowscalar(false)
        @test Array(logdensities(MeasureBase.Adapt.adapt(JLArray, pm), JLArray(flatview(X)))) ≈ logdensities(pm, X)
    end

    @testset "support of powers with array-variate bases" begin
        mv = MeasureBase.AsMeasure{typeof(MvNormal(zeros(2), I(2)))}(MvNormal(zeros(2), I(2)))
        Xf = randn(2, 3)
        @test insupport(mv^3, Xf) === true
        @test insupport(mv^3, nestedview(Xf)) === true
        @test logdensities(mv, Xf) ≈ [logpdf(mv.obj, Xf[:, j]) for j in 1:3]
        @test logdensityof(mv^3, Xf) ≈ sum(logpdf(mv.obj, Xf))
    end

    @testset "powers of array Diracs" begin
        D = Dirac([1.0, 2.0])
        μ = D^3
        x = [1.0 1.0 1.0; 2.0 2.0 2.0]
        @test MeasureBase.mspace_ndims(typeof(μ)) == 2
        @test logdensityof(μ, x) == 0 && insupport(μ, x)
        @test logdensityof(μ, [[1.0, 2.0] for _ in 1:3]) == 0
        @test logdensityof(μ, 2 .* x) == -Inf
        @test logdensities(μ, cat(x, 2 .* x; dims = 3)) == [0.0, -Inf]
        @test rand(μ) == nestedview(x)
    end

    @testset "streams with value-dependent sizes inside powers" begin
        bnd = mbind(a -> StdNormal()^(a > 0.5 ? 2 : 1), StdUniform(), vcat)
        @test testvalue(bnd) isa AbstractVector
        m = mcombine(vcat, StdNormal(), bnd^2)
        x = [0.5, 0.7, 0.1, 0.2, 0.3, 0.1]
        ℓ = logdensityof(m, x)
        @test ℓ ≈ logdensityof(StdNormal(), 0.5) + logdensityof(bnd, [0.7, 0.1, 0.2]) + logdensityof(bnd, [0.3, 0.1])
        @test logdensities(m, hcat(x, x)) ≈ [ℓ, ℓ]
        f = transport_to(StdUniform()^6, m)
        @test inverse(f)(f(x)) ≈ x
        @test flatview(inverse(f).(f.(sliced(hcat(x, x), Val(1))))) ≈ hcat(x, x)
    end

    @testset "powers of tuple products in streams" begin
        Pt = productmeasure((StdNormal(), StdExponential()^2))
        m = mcombine(vcat, StdNormal(), Pt^2)
        x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        ℓ = logdensityof(StdNormal(), 0.1) + logdensityof(Pt, (0.2, [0.3, 0.4])) + logdensityof(Pt, (0.5, [0.6, 0.7]))
        @test logdensityof(m, x) ≈ ℓ
        X = hcat(x, 2 .* x)
        @test logdensities(m, X) ≈ [logdensityof(m, X[:, j]) for j in 1:2]
        f = transport_to(StdUniform()^7, m)
        @test inverse(f)(f(x)) ≈ x
        @test flatview(inverse(f).(f.(sliced(X, Val(1))))) ≈ X
        Z, R = batched_transport_to_std_with_rest(StdUniform, Pt, hcat(x, x), (2,))
        @test size(Z) == (6, 2) && size(R) == (1, 2)
        Xr, _ = batched_transport_from_std_with_rest(StdUniform, Pt, Z, (2,))
        @test Xr[1] ≈ [0.1 0.1; 0.4 0.4] && size(Xr[2]) == (2, 2, 2)
    end

    @testset "batched kernels of measures without a declared rank" begin
        w = weightedmeasure(0.3, UnknownRankMeasure())
        @test logdensityof(w, [1.0, 2.0]) ≈ 0.3 - 3
        @test_throws ArgumentError MeasureBase.batched_logdensity_def(w, randn(2, 2))
        @test MeasureBase.batched_logdensity_def(weightedmeasure(0.3, StdNormal()^2), [1.0, 2.0]) === 0.3
    end

    @testset "ragged batches are evaluated variate by variate" begin
        V = VectorOfVectors([randn(2) for _ in 1:4])
        @test logdensities(StdNormal()^2, V) ≈ logdensityof.(Ref(StdNormal()^2), V)
        @test_throws ArgumentError logdensities(StdNormal(), VectorOfVectors([[1.0], [2.0], [3.0], [4.0]]))
        @test logdensities(StdNormal()^2, nestedview(randn(2, 4))) isa AbstractVector
    end

    @testset "static streams stay allocation-free" begin
        m1 = mcombine(vcat, StdNormal(), StdExponential()^static(2))
        x1 = SVector(0.1, 0.2, 0.3)
        @test logdensityof(m1, x1) ≈ logdensityof(StdNormal(), 0.1) + logdensityof(StdExponential()^2, [0.2, 0.3])
        @test @allocated(logdensityof(m1, x1)) == 0
        m2 = mcombine(vcat, StdNormal()^static(2), StdExponential()^static(3))
        x2 = SVector(0.1, 0.2, 0.3, 0.4, 0.5)
        @test @allocated(logdensityof(m2, x2)) == 0
        g = transport_to(StdUniform()^3, StdNormal()^3)
        v = randn(3)
        @test g.(SVector{3}(v))[] ≈ g(v)
        gs = transport_to(StdUniform(), StdNormal())
        @test gs.(SVector{3}(v)) isa SVector{3,Float64}
    end

    @testset "random variates of fused array products" begin
        P = productmeasure([pushfwd(Base.Fix1(*, s), StdExponential()) for s in 0.1:0.2:0.9])
        x = rand(GenContext{Float64}(), P)
        @test x isa Vector{Float64} && length(x) == 5
        X = batched_rand_impl(GenContext{Float64}(), P, (7,))
        @test size(X) == (5, 7)
        Pj = productmeasure(JLArray([weightedmeasure(log(i), StdNormal()) for i in 1:3]))
        @test length(rand(Pj)) == 3
    end
end
