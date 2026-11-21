# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# The samplers of discrete distributions for compute units other than the
# host: their kernels with given uniform variates, the distribution of
# their variates, and their use through `rand` on device arrays.
#
# The statistical tests may fail by chance with a probability of about
# 1e-8 (histograms) resp. 2e-9 (means and variances) each.

using Test
using Distributions, Random, Statistics
using MeasureBase
using MeasureBase: asmeasure, GenContext, batched_rand_impl, productmeasure, marginals
using MeasureBase.HeterogeneousComputing: AbstractComputeUnit, rootrng
import Adapt
using JLArrays
using StructArrays: StructArray

isdefined(@__MODULE__, :pmf_pvalue) || include("pmf_checks.jl")

@testset "samplers of discrete distributions" begin
    Ext = Base.get_extension(MeasureBase, :MeasureBaseDistributionsExt)
    JLArrays.allowscalar(false)
    jlunit = AbstractComputeUnit(JLBackend())
    jlctx(T, seed) = GenContext{T}(jlunit, rootrng(jlunit, seed))

    @testset "rejection driver" begin
        ctx = GenContext{Float64}(Xoshiro(1))
        # Round `i` proposes `i + param` and accepts the lanes of `script[i]`:
        script = [Bool[1, 0, 0, 0], Bool[0, 1, 0, 1], Bool[1, 1, 0, 0], Bool[0, 0, 1, 0]]
        for dev in (identity, JLArray)
            masks = map(dev, script)
            rounds = Ref(0)
            round = (_, param) -> (rounds[] += 1; (param .+ rounds[] .* one.(masks[rounds[]]), masks[rounds[]]))
            # Each lane keeps its first accepted candidate:
            x = MeasureBase._rejection_draw(round, ctx, true, 10)
            @test Array(x) == [11, 12, 14, 12] && rounds[] == 4
            # Inactive lanes take the first candidate and don't keep the loop going:
            rounds[] = 0
            x = MeasureBase._rejection_draw(round, ctx, dev(Bool[1, 1, 0, 1]), 10)
            @test Array(x) == [11, 12, 11, 12] && rounds[] == 2
        end
        # Single variates:
        rounds = Ref(0)
        @test MeasureBase._rejection_draw((_, p) -> (rounds[] += 1; (p + rounds[], rounds[] == 3)), ctx, true, 10) == 13
        # Constant RNGs repeat the first round, so it is the only one:
        rounds[] = 0
        cctx = GenContext{Float64}(MeasureBase.ConstantRNG())
        @test MeasureBase._rejection_draw((_, p) -> (rounds[] += 1; (p + rounds[], false)), cctx, true, 10) == 11
        @test Ext._poisson_draw(cctx, 200.7, (3,)) == fill(201, 3)
        @test Ext._poisson_draw(cctx, [3.2, 200.7], (2,)) == [3, 201]
    end

    @testset "Poisson kernels, $T" for T in (Float32, Float64)
        tol = 8 * eps(T) + 1e-13
        # Uniform variates in the unit interval, with both of its ends:
        grid = sort!(vcat(T[0, eps(T), prevfloat(one(T)), 1], rand(Xoshiro(2), T, 20_000)))

        # Inversion returns the quantile of the uniform variate, as far as
        # the precision resolves it:
        for λ in (0, 1e-3, 0.3, 1, 5, 9.99)
            d = Poisson(Float64(T(λ)))
            x = Ext._poisson_inversion.(T(λ), grid, eps(T) / 4)
            @test x isa Vector{Int} && issorted(x) && first(x) == 0
            @test all(cdf.(d, x .- 1) .<= grid .+ tol) && all(grid .<= cdf.(d, x) .+ tol)
        end

        kmax = Ext._max_count(T)
        for λ in (10, 14.048, 123.4, 1e5, Ext._poisson_max_rate(T))
            λT = T(λ)
            k = Ext._ptrs_candidate.(λT, grid, kmax)
            # No candidate at the ends of the unit interval and beyond the
            # exact integers, candidates increase with the variate:
            @test k isa Vector{Int} && first(k) == last(k) == -1 && all(>=(-1), k) && all(<=(kmax), k)
            @test issorted(filter(>=(0), k)) && count(>=(0), k) > 0.9 * length(k)
            @test Ext._ptrs_candidate(λT, one(T) / 2, kmax) == floor(Int, Float64(λT) - floor(Float64(λT)) + 0.43) + Int(floor(λT))
            # A zero second variate accepts every candidate, trials without
            # a candidate are never accepted:
            @test Ext._ptrs_accept.(λT, grid, zero(T), k) == (k .>= 0)
            @test !any(Ext._ptrs_accept.(λT, grid, one(T) / 2, -1))
        end

        # Candidates are exact integers also where the float type doesn't
        # resolve the sum of rate and offset:
        λT = Ext._poisson_max_rate(T)
        halfwidth = min(0.05, 2.0e5 / first(Ext._ptrs_constants(λT)))
        k = Ext._ptrs_candidate.(λT, T.(range(0.5 - halfwidth, 0.5 + halfwidth, length = 10^6)), kmax)
        @test sort!(unique(k)) == minimum(k):maximum(k)

        # Unsupported rates throw:
        ctx = GenContext{T}(Xoshiro(3))
        @test_throws "supports rates from 0 to" Ext._poisson_draw(ctx, 2 * kmax, (3,))
        @test_throws "supports rates from 0 to" Ext._poisson_draw(ctx, [1, -1], (2,))
        @test_throws "supports rates from 0 to" Ext._poisson_draw(ctx, [1, 2 * kmax], (2,))
    end

    @testset "PTRS envelope" begin
        # The acceptance probability of every trial is at most one, also at
        # rates where the envelope of the original algorithm falls short
        # (by 0.58 % at rate 14.048):
        us = range(0, 1, length = 400_001)[2:end-1]
        for λ in (10.0, 14.048, 30.86, 100.0, 1.0e3, 1.0e6)
            b, a, invα = Ext._ptrs_constants(λ)
            k = Ext._ptrs_candidate.(λ, us, Ext._max_count(Float64))
            hat = invα ./ (a ./ (0.5 .- abs.(us .- 0.5)) .^ 2 .+ b)
            @test maximum(pdf.(Poisson(λ), k) ./ hat) < 1
        end
    end

    @testset "distribution of Poisson variates, $T" for T in (Float32, Float64)
        n = 200_000
        ctx = GenContext{T}(Xoshiro(4))
        rates = [0.02, 0.9, 4.7, 9.99, 10, 12.3, 61.5, 1.0e3, 3.7e5, Ext._poisson_max_rate(T) / 2]
        for λ in rates
            d = Poisson(Float64(T(λ)))
            x = Ext._poisson_draw(ctx, λ, (n,))
            @test x isa Vector{Int}
            @test pmf_pvalue(x, d) > 1e-8
            @test abs(mean_zscore(x, d)) < 6 && abs(var_zscore(x, d)) < 6
        end
        # Arrays of rates in both regimes, and zero rates:
        mixed = [0, 0.9, 9.99, 10, 61.5, 1.0e3]
        X = Ext._poisson_draw(ctx, reshape(mixed, :, 1), (length(mixed), n))
        @test X isa Matrix{Int} && all(iszero, X[1, :])
        for (i, λ) in enumerate(mixed)
            @test pmf_pvalue(X[i, :], Poisson(Float64(T(λ)))) > 1e-8
        end
    end

    @testset "Poisson variates on device arrays" begin
        n = 100_000
        for T in (Float32, Float64), λ in (2.3, 47.5)
            m, d = asmeasure(Poisson(λ)), Poisson(Float64(T(λ)))
            x = batched_rand_impl(jlctx(T, 5), m, (n,))
            @test x isa JLArray{Int,1}
            @test pmf_pvalue(Array(x), d) > 1e-8
            @test abs(mean_zscore(Array(x), d)) < 6 && abs(var_zscore(Array(x), d)) < 6
            # The same stream key reproduces the draws, the RNG advances:
            ctx = jlctx(T, 5)
            @test Array(batched_rand_impl(ctx, m, (n,))) == Array(x)
            @test Array(batched_rand_impl(ctx, m, (n,))) != Array(x)
            @test rand(jlctx(T, 6), m) isa Int
            @test rand(jlctx(T, 6), m^5) isa JLArray{Int,1}
            @test size(batched_rand_impl(jlctx(T, 6), m^(2, 3), (4,))) == (2, 3, 4)
            @test size(batched_rand_impl(jlctx(T, 6), m, (0,))) == (0,)
        end

        # Array products draw with the rates of their marginals as a
        # parameter array, on the host or on the compute unit:
        rates = [0, 0.9, 9.99, 10, 61.5, 1.0e3]
        p_host = productmeasure(asmeasure.(Poisson.(rates)))
        p_dev = Adapt.adapt(jlunit, p_host)
        @test marginals(p_dev).obj.λ isa JLArray
        for p in (p_host, p_dev)
            ctx = jlctx(Float32, 7)
            x = rand(ctx, p)
            @test x isa JLArray{Int,1} && size(x) == (6,)
            X = batched_rand_impl(ctx, p, (n,))
            @test X isa JLArray{Int,2} && size(X) == (6, n)
            for (i, λ) in enumerate(rates)
                @test pmf_pvalue(Array(X)[i, :], Poisson(Float64(Float32(λ)))) > 1e-8
            end
        end
        p_mat = productmeasure(asmeasure.(Poisson.([1.0 20.0; 3.0 400.0])))
        @test size(rand(jlctx(Float32, 8), p_mat)) == (2, 2)
        @test size(batched_rand_impl(jlctx(Float32, 8), p_mat, (3, 2))) == (2, 2, 3, 2)

        # Host contexts draw the marginals with the samplers of Distributions:
        host_rates = [0.9, 20.0, 300.0]
        p = productmeasure(asmeasure.(Poisson.(host_rates)))
        draws(rng) = [rand(rng, Poisson(λ)) for λ in host_rates]
        @test rand(GenContext{Float32}(Xoshiro(10)), p) == draws(Xoshiro(10))
        @test batched_rand_impl(GenContext{Float32}(Xoshiro(10)), p, (4,)) ==
              (rng = Xoshiro(10); reduce(hcat, [draws(rng) for _ in 1:4]))

        # Rates beyond the exact integers of the context precision are refused:
        @test_throws "supports rates from 0 to" rand(jlctx(Float32, 9), asmeasure(Poisson(2.0e7)))
        @test_throws "supports rates from 0 to" rand(jlctx(Float32, 9), productmeasure(asmeasure.(Poisson.([1.0, 2.0e7]))))
        @test rand(jlctx(Float64, 9), asmeasure(Poisson(2.0e7))) isa Int
    end
end
