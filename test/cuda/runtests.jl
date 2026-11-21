# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# CUDA tests, not part of the default test suite (they need a CUDA GPU).
# From the repository root, Julia 1.11+ uses the local MeasureBase and the
# HeterogeneousComputing source revision declared in test/cuda/Project.toml:
#   julia --project=test/cuda -e 'using Pkg; Pkg.instantiate()'
#   julia --project=test/cuda test/cuda/runtests.jl

using Test
using CUDA
using Adapt: adapt
using HeterogeneousComputing: GenContext, AbstractComputeUnit, rootrng
using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, Dirac, asmeasure
using MeasureBase: productmeasure, pushfwd, mcombine, weightedmeasure, superpose
using MeasureBase: batched_rand_impl
using MeasureBase.InverseFunctions: inverse
using ArraysOfArrays: sliced, flatview
using AffineMaps: Mul, MulAdd
using Distributions: Normal, Uniform, Exponential, Logistic, Cauchy, Laplace, LogNormal, Weibull, Gamma, Beta
using Distributions: Poisson, Bernoulli, Binomial, MvNormal, Dirichlet
import PDMats

CUDA.allowscalar(false)

include(joinpath(@__DIR__, "..", "distributions", "pmf_checks.jl"))

# Evaluates `f` on device copies of `args` and compares with the plain
# result, the result must live on the device:
function test_cuda(f, args...)
    expected = f(args...)
    result = f(map(cu_copy, args)...)
    @test _device_array(result)
    @test _plain(result) ≈ _plain(expected) nans = true
    return result
end

cu_copy(x::AbstractArray) = CuArray(x)
cu_copy(μ::AbstractMeasure) = adapt(CuArray, μ)
cu_copy(x) = x
_device_array(x::AbstractArray) = parent_array(x) isa CuArray
_device_array(x::Tuple) = all(_device_array, x)
parent_array(x::CuArray) = x
parent_array(x::AbstractArray) = parent_array(parent(x))
parent_array(x::Base.ReshapedArray) = parent_array(parent(x))
_plain(x::AbstractArray) = Array(flatview(x))
_plain(x::Tuple) = map(_plain, x)

@testset "CUDA" begin
    X = randn(3, 20)
    Xc = vcat(randn(2, 20), rand(1, 20))

    @testset "densities" begin
        test_cuda(X -> logdensities(StdNormal()^3, X), X)
        test_cuda(X -> logdensities(StdNormal()^3, sliced(X, Val(1))), X)
        test_cuda(X -> logdensities((StdNormal()^3)^4, reshape(X[:, 1:16], 3, 4, 4)), X)
        test_cuda(X -> logdensities(weightedmeasure(log(0.3), StdNormal()^3), X), X)
        mc = mcombine(vcat, StdNormal()^2, StdUniform()^1)
        test_cuda(X -> logdensities(mc, X), Xc)
        mix = superpose(weightedmeasure(log(0.3), StdNormal()), weightedmeasure(log(0.7), StdUniform()))
        test_cuda(x -> logdensities(mix, x), rand(20))
        test_cuda(x -> logdensities(0.2 * StdNormal() + 0.8 * Dirac(0.0), x), vcat(randn(19), 0.0))
        test_cuda(x -> logdensities(Dirac(0.5), x), vcat(rand(19), 0.5))
        νe = pushfwd(Base.BroadcastFunction(exp), StdNormal()^3)
        test_cuda(Y -> logdensities(νe, Y), exp.(X))
        P = productmeasure([pushfwd(Mul(s), StdNormal()) for s in (1.0, 2.0, 3.0)])
        test_cuda((P, X) -> logdensities(P, X), P, X)
        test_cuda((P, X) -> logdensities(P, sliced(X, Val(1))), P, X)
        Pa = productmeasure([MeasureBase.AsMeasure{Normal{Float64}}(Normal(μ, 1.0)) for μ in (0.0, 1.0, 2.0)])
        test_cuda((P, X) -> logdensities(P, X), Pa, X)
        # Array products of array-variate marginals loop over the marginals
        # on the host, the marginals stay host arrays:
        Pv = productmeasure([weightedmeasure(log(i), StdNormal()^2) for i in 1:3])
        test_cuda(X -> logdensities(Pv, X), randn(2, 3, 5))

        q = mintegrate_exp(x -> -sum(abs2, x), StdNormal()^2)
        W = Float32[0.1 -0.5 1.0; 0.3 0.2 2.0]
        ℓ = test_cuda(X -> logdensities(q, X), W)
        @test Array(ℓ) ≈ -Float32(log(2π)) .- 1.5f0 .* vec(sum(abs2, W; dims = 1))

        r = restrict(x -> sum(x) > 0, StdNormal()^2)
        Xr = Float32[0.25 -1.5 0.5 -2.0; 0.5 0.25 -0.5 3.0]
        expected = -Float32(log(2π)) .- 0.5f0 .* vec(sum(abs2, Xr; dims = 1))
        expected[2:3] .= -Inf32
        ℓ = test_cuda(X -> logdensities(r, X), reshape(Xr, 2, 2, 2))
        @test size(ℓ) == (2, 2)
        @test vec(Array(ℓ)) ≈ expected

        mx = superpose(0.4f0 * StdUniform()^2, 0.6f0 * StdExponential()^2)
        Xm = Float32[0.5 -0.5 -1.0 2.5; 0.25 0.25 -0.5 0.25]
        ℓ = test_cuda(X -> logdensities(mx, X), Xm)
        @test Array(ℓ) ≈ Float32[log(0.4f0 + 0.6f0 * exp(-0.75f0)), -Inf, -Inf, log(0.6f0) - 2.75f0]
    end

    @testset "transports" begin
        g = transport_to(StdUniform()^3, StdNormal()^3)
        test_cuda(X -> flatview(g.(sliced(X, Val(1)))), X)
        test_cuda(X -> flatview(inverse(g).(g.(sliced(X, Val(1))))), X)
        h = transport_to(StdUniform()^(2, 3), (StdNormal()^2)^3)
        test_cuda(X -> flatview(h.(sliced(X, Val(2)))), randn(2, 3, 4))
        mc = mcombine(vcat, StdNormal()^2, StdUniform()^1)
        c = transport_to(StdExponential()^3, mc)
        test_cuda(X -> flatview(c.(sliced(X, Val(1)))), Xc)
        test_cuda(X -> flatview(inverse(c).(c.(sliced(X, Val(1))))), Xc)
        P = productmeasure([pushfwd(Mul(s), StdNormal()) for s in (1.0, 2.0, 3.0)])
        f = transport_to(StdUniform()^3, P)
        test_cuda((P, X) -> flatview(transport_to(StdUniform()^3, P).(sliced(X, Val(1)))), P, X)
        test_cuda((P, X) -> flatview(transport_to(P, StdUniform()^3).(sliced(X, Val(1)))), P, rand(3, 20))
        Pv = productmeasure([weightedmeasure(log(i), StdNormal()^2) for i in 1:3])
        fv = transport_to(StdUniform()^6, Pv)
        test_cuda(X -> flatview(fv.(sliced(X, Val(2)))), randn(2, 3, 5))
        test_cuda(X -> flatview(inverse(fv).(fv.(sliced(X, Val(2))))), randn(2, 3, 5))
        νe = pushfwd(Base.BroadcastFunction(exp), StdNormal()^3)
        fe = transport_to(StdUniform()^3, νe)
        test_cuda(Y -> flatview(fe.(sliced(Y, Val(1)))), exp.(X))
        A = [2.0 0.5; 0.0 1.5]
        b = [1.0, -1.0]
        νa = pushfwd(MulAdd(A, b), StdNormal()^2)
        νac = pushfwd(MulAdd(CuArray(A), CuArray(b)), StdNormal()^2)
        Ya = randn(2, 20)
        ra = flatview(transport_to(StdUniform()^2, νac).(sliced(CuArray(Ya), Val(1))))
        @test ra isa CuArray && Array(ra) ≈ flatview(transport_to(StdUniform()^2, νa).(sliced(Ya, Val(1))))
        # AffineMaps has no device support for log-abs-det-Jacobians yet:
        @test_broken Array(logdensities(νac, CuArray(Ya))) ≈ logdensities(νa, Ya)
    end

    @testset "wrapped distributions" begin
        for d in (Normal(0.3, 1.7), Uniform(-1.0, 2.5), Exponential(0.7), Logistic(0.2, 1.3), Cauchy(0.1, 0.8), Laplace(-0.4, 1.1), LogNormal(0.2, 0.6), Weibull(1.4, 0.9), Gamma(2.3, 1.2), Beta(2.5, 3.5))
            m = asmeasure(d)
            xd = rand(d, 20)
            test_cuda(x -> logdensities(m, x), xd)
            f = transport_to(StdNormal(), m)
            test_cuda(x -> f.(x), xd)
            test_cuda(z -> inverse(f).(z), randn(20))
        end
        for d in (Poisson(2.7), Bernoulli(0.3))
            test_cuda(x -> logdensities(asmeasure(d), x), Float64.(rand(d, 20)))
        end
        # Tails of the Float32 kernels:
        test_cuda(p -> transport_to(asmeasure(Cauchy(0f0, 1f0)), StdUniform()).(p), Float32[1f-30, 1f-9, 0.3, 0.5, 0.7, 1 - 6f-8])
        test_cuda(x -> logdensities(asmeasure(Logistic(0f0, 1f0)), x), Float32[-floatmax(Float32), -80, 0, 80, floatmax(Float32)])
        test_cuda(x -> logdensities(asmeasure(Cauchy(0f0, 1f0)), x), Float32[-1f30, -3, 0, 1f20, floatmax(Float32)])
        for Σ in ([1.7 0.5; 0.5 2.3], PDMats.PDiagMat([0.5, 2.0]), PDMats.ScalMat(2, 1.5))
            mvn = MvNormal([0.3, -2.9], Σ)
            mm = asmeasure(mvn)
            Xm = rand(mvn, 20)
            test_cuda((m, X) -> logdensities(m, X), mm, Xm)
            test_cuda((m, X) -> flatview(transport_to(StdNormal()^2, m).(sliced(X, Val(1)))), mm, Xm)
            test_cuda((m, Z) -> flatview(transport_to(m, StdNormal()^2).(sliced(Z, Val(1)))), mm, randn(2, 20))
        end
        # Single-precision kernels stay in single precision on the device:
        @test logdensities(asmeasure(Normal(0.3f0, 1.7f0)), CuArray(rand(Float32, 8))) isa CuArray{Float32}
        @test logdensities(cu_copy(asmeasure(MvNormal(Float32[0.3, -2.9], PDMats.PDiagMat(Float32[0.5, 2.0])))), CuArray(rand(Float32, 2, 8))) isa CuArray{Float32}
        dir = Dirichlet([2.0, 3.0, 4.0, 1.5])
        md = asmeasure(dir)
        Xd = rand(dir, 20)
        test_cuda((m, X) -> logdensities(m, X), md, Xd)
        test_cuda((m, X) -> flatview(transport_to(StdUniform()^3, m).(sliced(X, Val(1)))), md, Xd)
    end

    @testset "random variates" begin
        ctx = GenContext{Float32}(AbstractComputeUnit(CUDA.device()))
        for μ in (StdNormal(), StdUniform(), StdExponential(), StdNormal()^3, (StdNormal()^2)^3, weightedmeasure(0.3, StdNormal()^2), mcombine(vcat, StdNormal()^2, StdUniform()^1), pushfwd(Base.BroadcastFunction(exp), StdNormal()^2), superpose(StdNormal(), StdUniform()), 0.5 * StdNormal() + 0.5 * Dirac(0.0f0), asmeasure(Normal(0.3, 1.7)), asmeasure(Weibull(1.4, 0.9)), cu_copy(asmeasure(MvNormal([0.3, -2.9], [1.7 0.5; 0.5 2.3]))))
            X = batched_rand_impl(ctx, μ, (100,))
            @test X isa CuArray{Float32}
            ℓ = logdensities(μ, X)
            @test ℓ isa CuArray && all(isfinite, Array(ℓ))
        end
        Pt = productmeasure((StdNormal(), StdExponential()^2))
        Xt = batched_rand_impl(ctx, Pt, (50,))
        @test Xt isa Tuple && all(x -> x isa CuArray{Float32}, Xt)
        @test size(Xt[2]) == (2, 50)

        # Poisson variates are sampled on the device, by inversion for
        # small rates and by rejection for larger ones. The statistical
        # tests may fail by chance with a probability of about 1e-8 each:
        for λ in (2.7, 250.0)
            d = Poisson(Float64(Float32(λ)))
            x = batched_rand_impl(ctx, asmeasure(Poisson(λ)), (10^5,))
            @test x isa CuArray{Int,1}
            @test pmf_pvalue(Array(x), d) > 1e-8
            @test abs(mean_zscore(Array(x), d)) < 6 && abs(var_zscore(Array(x), d)) < 6
            @test rand(ctx, asmeasure(Poisson(λ))) isa Int
        end
        # Array products draw with the rates of their marginals as a
        # parameter array on the device:
        rates = [0, 0.9, 9.99, 10, 61.5, 1.0e3]
        Pp = cu_copy(productmeasure(asmeasure.(Poisson.(rates))))
        Xp = batched_rand_impl(ctx, Pp, (10^5,))
        @test Xp isa CuArray{Int,2} && size(Xp) == (6, 10^5)
        for (i, λ) in enumerate(rates)
            @test pmf_pvalue(Array(Xp)[i, :], Poisson(Float64(Float32(λ)))) > 1e-8
        end
        @test rand(ctx, Pp) isa CuArray{Int,1}
        @test_throws "supports rates from 0 to" rand(ctx, asmeasure(Poisson(2.0e7)))

        # Distributions without a standard transport or a sampler of their
        # own are refused on the device:
        @test_throws "host compute units" batched_rand_impl(ctx, asmeasure(Binomial(5, 0.3)), (100,))

        # RNGs with stream keys reproduce the draws:
        cunit = AbstractComputeUnit(CUDA.device())
        keyed_ctx(seed) = GenContext{Float32}(cunit, rootrng(cunit, seed))
        Xn = batched_rand_impl(keyed_ctx(1), StdNormal(), (100,))
        @test Array(batched_rand_impl(keyed_ctx(1), StdNormal(), (100,))) == Array(Xn)
        @test Array(batched_rand_impl(keyed_ctx(2), StdNormal(), (100,))) != Array(Xn)
    end
end
