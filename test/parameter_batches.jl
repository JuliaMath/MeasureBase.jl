# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Parameter batches: arrays of measures, likelihoods over arrays of
# parameter values and density measures over batches of variates, against
# the elementwise evaluation and independent formulas.

using Test
using MeasureBase
using MeasureBase: logdensities, likelihoodof, mintegrate, mintegrate_exp, productmeasure, pushfwd
using MeasureBase: weightedmeasure, restrict, asmeasure, StdNormal, StdExponential
using AffineMaps: MulAdd
import StructArrays
using StructArrays: StructArray
using ArraysOfArrays: VectorOfSimilarVectors
import JLArrays
using JLArrays: JLArray

@testset "parameter batches" begin
    normal_ld(x, μ, σ) = -((x - μ) / σ)^2 / 2 - log(σ) - log(2π) / 2
    normal(μ, σ) = pushfwd(MulAdd(σ, μ), StdNormal())

    μs, σs = [-0.3, 0.0, 0.4, 1.2, 2.0], [0.5, 1.0, 1.5, 0.8, 2.5]
    ms = normal.(μs, σs)
    ms_sa = StructArray(ms)

    @testset "arrays of measures with scalar variates" begin
        x, xs, X = 0.7, [0.1, -0.4, 0.9, 1.1, 3.0], [0.1 -0.4; 0.9 1.1; 3.0 0.2; -1.0 0.5; 2.2 0.0]
        for ms_ in (ms, ms_sa)
            @test logdensityof.(ms_, x) ≈ normal_ld.(x, μs, σs)
            @test logdensityof.(ms_, xs) ≈ normal_ld.(xs, μs, σs)
            @test logdensityof.(ms_, X) ≈ normal_ld.(X, μs, σs)
            @test logdensityof.(ms_, Ref(x)) ≈ normal_ld.(x, μs, σs)
            @test logdensityof.(ms_, permutedims(xs)) ≈ normal_ld.(permutedims(xs), μs, σs)
            @test densityof.(ms_, xs) ≈ exp.(normal_ld.(xs, μs, σs))
            @test logdensities(ms_, X) ≈ normal_ld.(X, μs, σs)
            @test logdensities(ms_, x) ≈ normal_ld.(x, μs, σs)
            # Fused with other broadcasts:
            @test logdensityof.(ms_, xs .+ 1) .- 1 ≈ normal_ld.(xs .+ 1, μs, σs) .- 1
        end
        # Support masks per measure:
        es = [weightedmeasure(w, StdExponential()) for w in (0.0, 0.5, 1.0)]
        @test logdensityof.(es, [1.0, -1.0, 2.0]) == [-1.0, -Inf, -1.0]
        @test logdensities(StructArray(es), [1.0 -1.0; -1.0 0.5; 2.0 3.0]) == [-1.0 -Inf; -Inf 0.0; -1.0 -2.0]
    end

    @testset "empty batches and single precision" begin
        @test logdensityof.(ms[1:0], 0.3) == Float64[]
        @test logdensities(ms_sa[1:0], zeros(0, 3)) == zeros(0, 3)
        ms32 = normal.(Float32.(μs), Float32.(σs))
        ℓ32 = logdensityof.(StructArray(ms32), 0.7f0)
        @test eltype(ℓ32) == Float32 && ℓ32 ≈ normal_ld.(0.7, μs, σs)
        @test eltype(logdensities(ms32, Float32[0.1, 0.2, 0.3, 0.4, 0.5])) == Float32
    end

    @testset "arrays of measures with array variates" begin
        n = 4
        pws = [m^n for m in ms]
        x = [0.2, -0.5, 1.3, 0.8]
        ref_shared = [sum(normal_ld.(x, μ, σ)) for (μ, σ) in zip(μs, σs)]
        @test logdensities(pws, x) ≈ ref_shared
        @test logdensities(StructArray(pws), x) ≈ ref_shared
        @test logdensityof.(pws, Ref(x)) ≈ ref_shared
        X = reshape(collect(range(-2, 2, length = n * 5 * 2)), n, 5, 2)
        ref_batch = [sum(normal_ld.(X[:, i, j], μs[i], σs[i])) for i in 1:5, j in 1:2]
        @test logdensities(pws, X) ≈ ref_batch
        @test logdensities(StructArray(pws), X) ≈ ref_batch
        @test logdensityof.(pws, VectorOfSimilarVectors(X[:, :, 1])) ≈ ref_batch[:, 1]
        @test_throws Exception logdensities(pws, X[1:3, :, :])
        # Powers of different shapes don't fuse, measure by measure:
        mixed = [ms[1]^2, ms[2]^3]
        @test logdensityof.(mixed, [[0.1, 0.2], [0.3, 0.4, 0.5]]) ≈
            [sum(normal_ld.([0.1, 0.2], μs[1], σs[1])), sum(normal_ld.([0.3, 0.4, 0.5], μs[2], σs[2]))]
    end

    @testset "likelihoods over parameter batches" begin
        data = [0.3, -1.2, 0.8, 2.1, 0.0, 1.5]
        Θ = StructArray((μ = μs, σ = σs))
        ref = [sum(normal_ld.(data, θ.μ, θ.σ)) for θ in Θ]
        ℒ = likelihoodof(θ -> normal(θ.μ, θ.σ)^length(data), data)
        @test logdensityof.(Ref(ℒ), Θ) ≈ ref
        @test logdensityof.(ℒ, Θ) ≈ ref
        @test logdensityof.(Ref(ℒ), Θ) ≈ [logdensityof(ℒ, θ) for θ in Θ]
        @test densityof.(Ref(ℒ), Θ) ≈ [densityof(ℒ, θ) for θ in Θ]
        @test logdensityof.(Ref(ℒ), Θ[1:0]) == Float64[]
        # Scalar observations and plain arrays of parameters:
        ℒ1 = likelihoodof(σ -> normal(0.5, σ), 1.2)
        @test logdensityof.(Ref(ℒ1), σs) ≈ normal_ld.(1.2, 0.5, σs)
        # Structured observations are evaluated measure by measure:
        ℒt = likelihoodof(θ -> productmeasure((a = normal(θ.μ, 1.0), b = normal(0.0, θ.σ))), (a = 0.4, b = -0.3))
        @test logdensityof.(Ref(ℒt), Θ) ≈ normal_ld.(0.4, μs, 1.0) .+ normal_ld.(-0.3, 0.0, σs)
    end

    @testset "density measures over batches of variates" begin
        data = [0.3, -1.2, 0.8]
        ℒ = likelihoodof(θ -> normal(θ.μ, θ.σ)^length(data), data)
        prior = productmeasure((μ = StdNormal(), σ = StdExponential()))
        post = mintegrate(ℒ, prior)
        Θ = StructArray((μ = [0.1, -0.5, 1.0, 0.2], σ = [0.5, 1.5, 2.0, -1.0]))
        ref = [logdensityof(post, θ) for θ in Θ]
        @test ref[4] == -Inf
        @test logdensities(post, Θ) ≈ ref
        @test logdensities(post, (μ = Θ.μ, σ = Θ.σ)) ≈ ref
        @test isfinite(ref[1]) && ref[1] ≈ sum(normal_ld.(data, 0.1, 0.5)) + normal_ld(0.1, 0, 1) - 0.5
        # The integrand is not evaluated outside of the support of the base:
        guarded = likelihoodof(θ -> (θ.σ > 0 || error("σ ≤ 0")) ? normal(θ.μ, θ.σ) : nothing, 0.3)
        @test logdensities(mintegrate(guarded, prior), Θ)[4] == -Inf
        # Vector variates and log-density integrands:
        post_v = mintegrate_exp(x -> -sum(abs2, x), StdNormal()^2)
        Xv = [0.1 0.5 -1.0; 0.3 -0.2 2.0]
        @test logdensities(post_v, Xv) ≈ [logdensityof(post_v, Xv[:, j]) for j in 1:3]
        @test logdensities(post_v, Xv) ≈ [-1.5 * sum(abs2, Xv[:, j]) - log(2π) for j in 1:3]
        post_w = mintegrate(ℒ, weightedmeasure(0.5, prior))
        @test logdensities(post_w, Θ) ≈ [logdensityof(post_w, θ) for θ in Θ]
        @test logdensities(post, Θ[1:0]) == Float64[]
        # Vector variates outside of the support of the base:
        post_e = mintegrate_exp(x -> -sum(x), StdExponential()^2)
        Xe = [0.1 -0.5 1.0; 0.3 0.2 2.0]
        @test logdensities(post_e, Xe) ≈ [-0.8, -Inf, -6.0]
    end

    @testset "device arrays" begin
        JLArrays.allowscalar(false)
        ms_j = StructArray(ms)
        ms_j = StructArray{eltype(ms)}(map(JLArray, StructArrays.components(ms_j)))
        xs = [0.1, -0.4, 0.9, 1.1, 3.0]
        ℓ = logdensityof.(ms_j, JLArray(xs))
        @test ℓ isa JLArray && Array(ℓ) ≈ normal_ld.(xs, μs, σs)
        data = [0.3, -1.2, 0.8, 2.1]
        ℒ = likelihoodof(θ -> normal(θ.μ, θ.σ)^4, JLArray(data))
        Θ = StructArray((μ = JLArray(μs), σ = JLArray(σs)))
        ℓθ = logdensityof.(Ref(ℒ), Θ)
        @test ℓθ isa JLArray && Array(ℓθ) ≈ [sum(normal_ld.(data, μ, σ)) for (μ, σ) in zip(μs, σs)]
        post = mintegrate_exp(x -> -x, StdExponential())
        ℓp = logdensities(post, JLArray([0.5, -1.0, 2.0]))
        @test ℓp isa JLArray && Array(ℓp) == [-1.0, -Inf, -4.0]
    end
end
