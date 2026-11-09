# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random, Statistics
using StableRNGs: StableRNG
using StructArrays: StructArray
using Adapt: adapt
using JLArrays
using AffineMaps: Mul
using ArraysOfArrays: sliced, flatview
using InverseFunctions: inverse

using MeasureBase
using MeasureBase: StdNormal, StdUniform, productmeasure, pushfwd, marginals, transport_to, logdensities
import Distributions

@testset "products over arrays of marginals" begin
    stblrng() = StableRNG(789990641)

    @testset "struct array storage" begin
        P = productmeasure([pushfwd(Mul(s), StdNormal()) for s in (1.0, 2.0, 3.0)])
        mar = marginals(P)
        @test mar isa StructArray
        @test length(mar) == 3 && mar[2] == pushfwd(Mul(2.0), StdNormal())
        @test productmeasure(mar) == P

        x = randn(stblrng(), 3)
        ℓ_ref = sum(logdensityof(m, xi) for (m, xi) in zip(mar, x))
        @test @inferred(logdensityof(P, x)) ≈ ℓ_ref
        @test @inferred(MeasureBase.logdensity_def(P, x)) isa Real
        X = randn(stblrng(), 3, 5)
        @test @inferred(logdensities(P, X)) ≈ [logdensityof(P, x) for x in eachcol(X)]
        @test logdensities(P, sliced(X, Val(1))) ≈ logdensities(P, X)

        f = transport_to(StdUniform()^3, P)
        @test inverse(f)(f(x)) ≈ x
        @test flatview(f.(sliced(X, Val(1)))) ≈ stack(map(f, eachcol(X)))
        @test flatview(inverse(f).(f.(sliced(X, Val(1))))) ≈ X

        @test rand(stblrng(), P) isa Vector{Float64}
        Xr = rand(stblrng(), P^100)
        @test size(flatview(Xr)) == (3, 100)
        @test isapprox(vec(mean(flatview(Xr), dims = 2)), zeros(3), atol = 1.0)
        @test isapprox(vec(std(flatview(Xr), dims = 2)), [1.0, 2.0, 3.0], rtol = 0.4)

        # Non-isbits marginals keep their container:
        Pv = productmeasure([pushfwd(Mul(randn(stblrng(), 2, 2)), StdNormal()^2) for _ in 1:2])
        @test !(marginals(Pv) isa StructArray)
    end

    @testset "marginals without a constructor from their fields" begin
        # InverseGamma stores a Gamma but has no constructor taking one:
        d = Distributions.InverseGamma(2.0, 1.0)
        ld(x) = sum(xi -> Distributions.logpdf(d, xi), x)
        for P in (productmeasure([asmeasure(d) for _ in 1:2]), asmeasure(Distributions.product_distribution(fill(d, 2))))
            @test logdensityof(P, [1.0, 2.0]) ≈ ld([1.0, 2.0])
            @test logdensities(P, [1.0 0.5; 2.0 3.0]) ≈ [ld([1.0, 2.0]), ld([0.5, 3.0])]
        end
    end

    @testset "draws over mapped arrays are eager" begin
        mapped = MeasureBase.MappedArrays.mappedarray(n -> StdNormal()^n, [1, 2])
        x = rand(stblrng(), MeasureBase.ProductMeasure(mapped))
        @test x[1] == x[1] && x[2] == x[2]
        @test length.(x) == [1, 2]
    end

    @testset "device arrays" begin
        JLArrays.allowscalar(false)
        P = productmeasure([pushfwd(Mul(s), StdNormal()) for s in (1.0, 2.0, 3.0)])
        Pj = adapt(JLArray, P)
        @test marginals(Pj) isa StructArray
        X = randn(stblrng(), 3, 5)
        Xj = JLArray(X)
        ℓj = logdensities(Pj, Xj)
        @test ℓj isa JLArray && Array(ℓj) ≈ logdensities(P, X)
        f = transport_to(StdUniform()^3, Pj)
        Yj = f.(sliced(Xj, Val(1)))
        @test flatview(Yj) isa JLArray
        @test Array(flatview(Yj)) ≈ flatview(transport_to(StdUniform()^3, P).(sliced(X, Val(1))))
    end
end
