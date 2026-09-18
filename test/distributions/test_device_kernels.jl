# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# The device-friendly kernels of the wrapped distribution families: their
# densities and transports agree with Distributions, and they run on flat
# batches, also of device arrays.

using Test
using Distributions, LinearAlgebra, StableRNGs, Statistics
using MeasureBase
using MeasureBase: asmeasure, GenContext, StdUniform, StdNormal, StdExponential, StdLogistic
using MeasureBase: batched_rand_impl, logdensities, insupport
using MeasureBase.InverseFunctions: inverse
using ArraysOfArrays: sliced, flatview
import Adapt
using JLArrays

@testset "device kernels of distribution families" begin
    JLArrays.allowscalar(false)
    stblrng() = StableRNG(28734)

    families = [
        Normal(0.3, 1.7), Uniform(-1.0, 2.5), Exponential(0.7), Logistic(0.2, 1.3), Cauchy(0.1, 0.8),
        Laplace(-0.4, 1.1), LogNormal(0.2, 0.6), Weibull(1.4, 0.9), Gamma(2.3, 1.2), Beta(2.5, 3.5),
        Normal(0.3f0, 1.7f0), Gamma(0.7, 2.0), Beta(0.6, 0.8),
    ]

    @testset "$(nameof(typeof(d)))" for d in families
        m = asmeasure(d)
        xs = vcat(rand(stblrng(), d, 20), [-1.0, 0.0, 1.0, 5.0, Inf])
        xs = eltype(d) == Float32 ? Float32.(xs) : xs
        ℓ_ref = logpdf.(d, xs)
        @test all(map((a, b) -> a == b || a ≈ b || (isnan(a) && isnan(b)), logdensityof.(Ref(m), xs), ℓ_ref))
        @test logdensities(m, xs) ≈ ℓ_ref nans = true
        @test Array(logdensities(m, JLArray(xs))) ≈ ℓ_ref nans = true
        @test insupport.(Ref(m), xs) == Distributions.insupport.(d, xs)

        if d isa ContinuousUnivariateDistribution
            x = rand(stblrng(), d, 12)
            f = transport_to(StdUniform(), m)
            p = f.(x)
            @test p ≈ cdf.(d, x)
            @test inverse(f).(p) ≈ x
            @test Array(inverse(f).(f.(JLArray(x)))) ≈ x
            g = transport_to(StdNormal(), m)
            @test inverse(g).(g.(x)) ≈ x
            X = batched_rand_impl(GenContext{Float64}(stblrng()), m, (5,))
            @test X isa Vector{Float64} && length(X) == 5
            # The transport-based generation used on devices, in single precision:
            Xs = MeasureBase._rand_default(GenContext{Float32}(stblrng()), m, (2000,), MeasureBase._NoRandImpl())
            @test Xs isa Vector{Float32} && all(insupport.(Ref(m), Xs))
            if isfinite(mean(d)) && isfinite(var(d))
                @test isapprox(mean(Xs), mean(d), atol = 5 * sqrt(var(d) / 2000) + 1e-3)
            end
        end
    end

    @testset "discrete families" begin
        for d in (Poisson(2.7), Bernoulli(0.3))
            m = asmeasure(d)
            xs = [0, 1, 2, 3, 7]
            @test logdensityof.(Ref(m), xs) ≈ logpdf.(d, xs)
            @test logdensities(m, xs) ≈ logpdf.(d, xs)
            @test Array(logdensities(m, JLArray(xs))) ≈ logpdf.(d, xs)
            @test logdensityof(m, -1) == -Inf && logdensityof(m, 1.5) == -Inf
        end
    end

    @testset "MvNormal" begin
        for Σ in [[1.7 0.5; 0.5 2.3], PDMats.PDiagMat([0.5, 2.0]), PDMats.ScalMat(2, 1.5)]
            d = MvNormal([0.3, -2.9], Σ)
            m = asmeasure(d)
            X = rand(stblrng(), d, 6)
            ℓ_ref = logpdf(d, X)
            @test logdensityof(m, X[:, 1]) ≈ ℓ_ref[1]
            @test logdensities(m, X) ≈ ℓ_ref
            @test logdensities(m, sliced(X, Val(1))) ≈ ℓ_ref
            @test logdensityof(m^6, X) ≈ sum(ℓ_ref)
            f = transport_to(StdNormal()^2, m)
            Y = f.(sliced(X, Val(1)))
            @test flatview(Y) ≈ stack(map(f, eachcol(X)))
            @test flatview(inverse(f).(Y)) ≈ X
            # JLArrays have no triangular solves, so only diagonal
            # covariances run on them (CUDA covers the general case):
            if !(Σ isa AbstractMatrix)
                mj = Adapt.adapt(JLArray, m)
                @test Array(logdensities(mj, JLArray(X))) ≈ ℓ_ref
                fj = transport_to(StdNormal()^2, mj)
                @test Array(flatview(fj.(sliced(JLArray(X), Val(1))))) ≈ flatview(Y)
                @test Array(flatview(inverse(fj).(sliced(JLArray(flatview(Y)), Val(1))))) ≈ X
            end
            @test size(batched_rand_impl(GenContext{Float64}(stblrng()), m, (7,))) == (2, 7)
        end
    end

    @testset "Dirichlet" begin
        d = Dirichlet([2.0, 3.0, 4.0, 1.5])
        m = asmeasure(d)
        X = rand(stblrng(), d, 6)
        ℓ_ref = logpdf(d, X)
        @test logdensityof(m, X[:, 1]) ≈ ℓ_ref[1]
        @test logdensities(m, X) ≈ ℓ_ref
        @test logdensityof(m, [0.5, 0.5, 0.2, 0.1]) == -Inf
        @test logdensityof(m, [0.5, 0.6, -0.1, 0.0]) == -Inf
        mj = Adapt.adapt(JLArray, m)
        @test Array(logdensities(mj, JLArray(X))) ≈ ℓ_ref
        f = transport_to(StdUniform()^3, m)
        Y = f.(sliced(X, Val(1)))
        @test flatview(Y) ≈ stack(map(f, eachcol(X)))
        @test flatview(inverse(f).(Y)) ≈ X
        @test all(0 .<= flatview(Y) .<= 1)
        Xr = batched_rand_impl(GenContext{Float64}(stblrng()), m, (200,))
        @test size(Xr) == (4, 200) && all(sum(Xr; dims = 1) .≈ 1)
    end
end
