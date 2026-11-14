# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# The device-friendly kernels of the wrapped distribution families: their
# densities and transports agree with Distributions, and they run on flat
# batches, also of device arrays.

using Test
using Distributions, LinearAlgebra, StableRNGs, Statistics
using MeasureBase
using MeasureBase: asmeasure, GenContext, StdUniform, StdNormal, StdExponential, StdLogistic
using MeasureBase: batched_rand_impl, logdensities, insupportof, static, productmeasure
using MeasureBase.HeterogeneousComputing: rootrng
using MeasureBase.InverseFunctions: inverse
using MeasureBase.ChangesOfVariables: with_logabsdet_jacobian
using ArraysOfArrays: sliced, flatview
import Adapt
using JLArrays
using StructArrays: StructArray

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
        @test insupportof.(Ref(m), xs) == (Distributions.insupport.(d, xs) .& isfinite.(xs))

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
            @test Xs isa Vector{Float32} && all(insupportof.(Ref(m), Xs))
            if isfinite(mean(d)) && isfinite(var(d))
                @test isapprox(mean(Xs), mean(d), atol = 5 * sqrt(var(d) / 2000) + 1e-3)
            end
        end
    end

    @testset "single precision kernels stay in single precision" begin
        for d in (Normal(0.3f0, 1.7f0), LogNormal(0.2f0, 0.6f0), Gamma(2.3f0, 1.2f0), Beta(2.5f0, 3.5f0), Poisson(2.7f0))
            m = asmeasure(d)
            xs = d isa Poisson ? [0, 1, 3, 4] : Float32[0.1, 0.4, 0.5, 0.9]
            @test logdensityof(m, xs[2]) isa Float32
            @test eltype(logdensities(m, xs)) == Float32
            @test eltype(logdensities(m, JLArray(xs))) == Float32
            @test eltype(logdensities(m^2, reshape(xs, 2, 2))) == Float32
        end
        for Σ in (PDMats.PDMat(Float32[1.7 0.5; 0.5 2.3]), PDMats.PDiagMat(Float32[0.5, 2.0]), PDMats.ScalMat(2, 1.5f0))
            m = asmeasure(MvNormal(Float32[0.3, -2.9], Σ))
            @test eltype(logdensities(m, rand(stblrng(), Float32, 2, 3))) == Float32
            @test eltype(logdensities(m, zeros(Float32, 2, 0))) == Float32
        end
        f = transport_to(StdNormal(), asmeasure(Normal(0.3f0, 1.7f0)))
        @test with_logabsdet_jacobian(f, 1f0) isa Tuple{Float32,Float32}
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

    @testset "discrete distributions without standard transports" begin
        jl32 = GenContext{Float32}(MeasureBase.HeterogeneousComputing.AbstractComputeUnit(JLBackend()))
        # Bernoulli variates compare uniform variates with the success
        # probability:
        X = batched_rand_impl(jl32, asmeasure(Bernoulli(0.3)), (10^4,))
        @test X isa JLArray{Bool,1}
        @test count(Array(X)) / 10^4 ≈ 0.3 atol = 5 * sqrt(0.21 / 10^4)
        @test rand(jl32, asmeasure(Bernoulli(0.3))) isa Bool
        @test rand(jl32, asmeasure(Bernoulli(0.3))^3) isa JLArray{Bool,1}
        @test all(Array(batched_rand_impl(jl32, asmeasure(Bernoulli(1.0)), (100,))))
        @test !any(Array(batched_rand_impl(jl32, asmeasure(Bernoulli(0.0)), (100,))))
    end

    @testset "host draws for other compute units" begin
        jlunit = MeasureBase.HeterogeneousComputing.AbstractComputeUnit(JLBackend())
        jl32 = GenContext{Float32}(jlunit)
        # Distributions without standard transports are sampled on the host
        # and moved to the compute unit, integer variates stay exact:
        for d in (Poisson(2.0), Categorical([0.2, 0.8]), Binomial(5, 0.3))
            m = asmeasure(d)
            @test rand(jl32, m) isa Int
            X = batched_rand_impl(jl32, m, (10^4,))
            @test X isa JLArray{Int,1}
            @test mean(Array(X)) ≈ mean(d) atol = 5 * std(d) / 100
            @test rand(jl32, m^static(3)) isa JLArray{Int,1}
            @test rand(jl32, m^(2, 3)) isa JLArray{Int,2}
        end
        mn = asmeasure(Multinomial(5, [0.2, 0.3, 0.5]))
        @test rand(jl32, mn) isa JLArray{Int,1}
        Xmn = flatview(rand(jl32, mn^4))
        @test Xmn isa JLArray{Int,2} && size(Xmn) == (3, 4) && all(==(5), sum(Array(Xmn); dims = 1))
        w = asmeasure(Wishart(3.0, [1.0 0.2; 0.2 1.0]))
        @test rand(jl32, w) isa JLArray{Float32,2}
        Xw = flatview(rand(jl32, w^4))
        @test Xw isa JLArray{Float32,3} && size(Xw) == (2, 2, 4)

        # The draws follow from the RNG of the context, which advances:
        keyed_ctx(seed) = GenContext{Float32}(jlunit, rootrng(jlunit, seed))
        pois = asmeasure(Poisson(20.0))^8
        ctx = keyed_ctx(1)
        x1, x2 = Array(rand(ctx, pois)), Array(rand(ctx, pois))
        @test x1 != x2
        ctx = keyed_ctx(1)
        @test Array(rand(ctx, pois)) == x1 && Array(rand(ctx, pois)) == x2
        @test Array(rand(keyed_ctx(2), pois)) != x1
        # Host contexts draw from their RNG directly:
        @test rand(GenContext{Float64}(stblrng()), pois) == rand(stblrng(), Poisson(20.0), 8)

        # Array products of such distributions draw on the host as a whole,
        # also with parameter columns on the compute unit:
        rates = Float32[0.5, 2.0, 30.0]
        counts = productmeasure(asmeasure.(Poisson.(rates)))
        for p in (counts, Adapt.adapt(jlunit, counts))
            @test rand(jl32, p) isa JLArray{Int,1}
            X = flatview(rand(keyed_ctx(3), p^(10^4)))
            @test X isa JLArray{Int,2} && size(X) == (3, 10^4)
            @test all(abs.(vec(mean(Array(X); dims = 2)) .- rates) .< 5 .* sqrt.(rates ./ 10^4))
        end
        @test rand(jl32, productmeasure(asmeasure.(Bernoulli.(Float32[0.1, 0.9])))) isa JLArray{Bool,1}
        # Products of array variates are arrays of variates on the compute
        # unit, flat batches of them are device arrays:
        pmn = productmeasure([mn, mn])
        @test rand(jl32, pmn) isa Vector{<:JLArray{Int,1}}
        @test batched_rand_impl(jl32, pmn, (4,)) isa JLArray{Int}
        # Products of distributions with standard transports keep drawing
        # on the compute unit, also with parameters on it:
        mvn = Adapt.adapt(jlunit, asmeasure(MvNormal(zeros(Float32, 2), Matrix{Float32}(I, 2, 2))))
        @test rand(jl32, productmeasure([mvn, mvn])) isa Vector{<:JLArray{Float32,1}}

        # Variates that can't be moved to the compute unit are refused:
        @test_throws "host compute units" rand(jl32, asmeasure(LKJCholesky(3, 2.0)))
        @test_throws "host compute units" rand(jl32, productmeasure(fill(asmeasure(LKJCholesky(3, 2.0)), 2)))
    end

    @testset "parameter batches of distributions" begin
        μs, σs = [0.3, -1.0, 2.0], [1.7, 0.5, 1.2]
        ds = Normal.(μs, σs)
        xs = [0.1, -0.4, 2.5]
        @test logdensityof.(asmeasure.(ds), xs) ≈ logpdf.(ds, xs)
        @test Array(logdensityof.(JLArray(asmeasure.(ds)), JLArray(xs))) ≈ logpdf.(ds, xs)
        data = [0.2, 0.9, -0.3]
        ℒ = likelihoodof(θ -> asmeasure(Normal(θ.μ, θ.σ))^3, JLArray(data))
        Θ = StructArray((μ = JLArray(μs), σ = JLArray(σs)))
        @test Array(logdensityof.(Ref(ℒ), Θ)) ≈ [sum(logpdf.(d, data)) for d in ds]
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
            # Triangular matrix products give wrong values on JLArrays, so
            # only diagonal and scalar covariances run on them (CUDA covers
            # the general case):
            if Σ isa Union{PDMats.PDiagMat,PDMats.ScalMat}
                mj = Adapt.adapt(JLArray, m)
                @test Array(logdensities(mj, JLArray(X))) ≈ ℓ_ref
                fj = transport_to(StdNormal()^2, mj)
                @test Array(flatview(fj.(sliced(JLArray(X), Val(1))))) ≈ flatview(Y)
                @test Array(flatview(inverse(fj).(sliced(JLArray(flatview(Y)), Val(1))))) ≈ X
                Xj = batched_rand_impl(GenContext{Float32}(MeasureBase.HeterogeneousComputing.AbstractComputeUnit(JLBackend())), mj, (7,))
                @test Xj isa JLArray && size(Xj) == (2, 7) && all(isfinite, Array(logdensities(mj, Xj)))
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

@testset "wrapped distributions outside the support" begin
    for d in (Uniform(-1.0, 2.5), Exponential(0.7), LogNormal(0.2, 0.6), Weibull(1.4, 0.9), Gamma(2.3, 1.2), Beta(2.5, 3.5))
        m = asmeasure(d)
        x = minimum(d) - 1
        @test logdensityof(m, x) == -Inf
        @test isnan(transport_to(StdNormal(), m)(x))
        @test isnan(transport_to(StdNormal(), m)(Inf))
        Y = Array(transport_to(StdNormal(), m).(JLArray([x, mean(d)])))
        @test isnan(Y[1]) && !isnan(Y[2])
    end
    for d in (Cauchy(0.1, 0.8), Laplace(-0.4, 1.1), Gamma(2.3, 1.2), Beta(2.5, 3.5), Uniform(-1.0, 2.5))
        f = transport_to(asmeasure(d), StdUniform())
        @test isnan(f(1.5)) && isnan(f(-0.5)) && !isnan(f(0.3))
        @test isequal(Array(f.(JLArray([1.5, 0.3]))), f.([1.5, 0.3]))
    end
    for d in (Exponential(0.7), Weibull(1.4, 0.9))
        @test isnan(transport_to(asmeasure(d), StdExponential())(-1.0))
    end
    for d in (Normal(0.3, 1.7), Logistic(0.2, 1.3), 2 * TDist(3.0) + 1), S in (StdNormal(), StdUniform(), StdLogistic()), x in (Inf, -Inf)
        @test isnan(transport_to(S, asmeasure(d))(x))
    end
    mn = asmeasure(MvNormal([0.0, 0.0], [1 0.5; 0.5 1]))
    for S in (StdNormal()^2, StdUniform()^2)
        f = transport_to(S, mn)
        @test all(isnan, f([Inf, 0.0])) && all(isnan, f([0.0, -Inf])) && !any(isnan, f([0.0, 0.0]))
        Z = flatview(f.(sliced([Inf 0.3; 0.0 0.1], Val(1))))
        @test all(isnan, Z[:, 1]) && !any(isnan, Z[:, 2])
    end
    @test logdensityof(asmeasure(Rayleigh(2.0)), Inf) == -Inf
    @test logdensityof(asmeasure(Poisson(2.7)), Inf) == -Inf
    @test logdensityof(asmeasure(Poisson(2.7)), 1.5) == -Inf
    @test logdensityof(asmeasure(Poisson(2.7)), -1.0) == -Inf
    @test logdensityof(asmeasure(Bernoulli(0.3)), 0.5) == -Inf
    md = asmeasure(Dirichlet([2.0, 3.0, 4.0, 1.5]))
    @test all(isnan, transport_to(StdUniform()^3, md)([0.5, 0.5, 0.5, 0.5]))
    @test all(isnan, transport_to(md, StdUniform()^3)([1.5, 0.5, 0.5]))
    @test logdensityof(md, [0.5, 0.5, 0.5, 0.5]) == -Inf
    Zd = flatview(transport_to(StdUniform()^3, md).(sliced([0.5 0.1; 0.5 0.2; 0.5 0.3; 0.5 0.4], Val(1))))
    @test all(isnan, Zd[:, 1]) && !any(isnan, Zd[:, 2])
end
