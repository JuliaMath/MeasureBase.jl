# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseTestExt

using Test

using MeasureBase
using MeasureBase: basemeasure_depth, proxy, insupportof, transport_to
import MeasureBase.Interface: test_interface, test_transport, test_smf

using DensityInterface: logdensityof
using InverseFunctions: inverse
using ChangesOfVariables: with_logabsdet_jacobian
using Tricks: static_hasmethod
using IntervalSets: Interval

function dynamic_basemeasure_depth(μ::M) where {M}
    if static_hasmethod(proxy, Tuple{M})
        π = proxy(μ)
        if static_hasmethod(basemeasure, Tuple{typeof(π)})
            basemeasure(π) == basemeasure(μ) && return dynamic_basemeasure_depth(π)
        end
    end
    β = basemeasure(μ)
    depth = 0
    while μ ≠ β
        (μ, β) = (β, basemeasure(β))
        depth += 1
    end
    return depth
end

function test_interface(μ::M) where {M}
    @eval begin
        μ = $μ
        @testset "$μ" begin
            μ = $μ

            @test AbstractMeasure(μ) isa AbstractMeasure

            ###########################################################################
            # basemeasure_depth
            static_depth = @inferred basemeasure_depth(μ)

            dynamic_depth = dynamic_basemeasure_depth(μ)

            if static_depth > dynamic_depth
                @warn "basemeasure_depth($μ) greater than requirement, could add some overhead"
            end
            @test static_depth ≥ dynamic_depth

            ###########################################################################
            # testvalue, logdensityof

            x = @inferred testvalue(Float64, μ)
            β = @inferred basemeasure(μ)

            ℓμ = @inferred logdensityof(μ, x)
            ℓβ = @inferred logdensityof(β, x)

            @test ℓμ ≈ logdensity_def(μ, x) + ℓβ

            @test logdensity_def(μ, testvalue(Float64, μ)) isa Real

            test_batched_densities(asmeasure(μ))
        end
    end
end

structisapprox(a, b) = isapprox(a, b)
structisapprox(a::NTuple{N,Any}, b::NTuple{N,Any}) where {N} = all(map(structisapprox, a, b))
function structisapprox(a::NamedTuple{names}, b::NamedTuple{names}) where {names}
    all(map(structisapprox, values(a), values(b)))
end

batchisapprox(A, B) = size(A) == size(B) && all(map(structisapprox, A, B))

# A batch of variates of `μ`: random draws, or copies of the test value if
# `μ` has infinite mass and so no normalized law to draw from.
test_batch(μ) = _test_batch(μ, massof(μ))
_test_batch(μ, mass::Real) = mass == Inf ? fill(testvalue(μ), 3) : rand(μ^3)
_test_batch(μ, ::Any) = rand(μ^3)

function test_batched_densities(μ)
    X = test_batch(μ)
    ℓ = map(x -> logdensityof(μ, x), X)
    @test size(X) == (3,) && all(>(-Inf), ℓ)
    @test batchisapprox(logdensities(μ, X), ℓ)
    @test logdensityof(μ^3, X) ≈ sum(ℓ)
    @test isempty(logdensities(μ, X[1:0]))
end

function test_transport(ν, μ)
    supertype(x) = Any
    supertype(x::Real) = Real
    supertype(x::AbstractArray{<:Real,N}) where {N} = AbstractArray{<:Real,N}

    @testset "transport_to $μ to $ν" begin
        x = rand(μ)
        f = transport_to(ν, μ)
        y = @inferred f(x)
        @test structisapprox(@inferred(inverse(f)(y)), x)
        @test @inferred(with_logabsdet_jacobian(f, x)) isa Tuple{supertype(y),Real}
        @test @inferred(with_logabsdet_jacobian(inverse(f), y)) isa Tuple{supertype(x),Real}
        y2, ladj_fwd = with_logabsdet_jacobian(f, x)
        x2, ladj_inv = with_logabsdet_jacobian(inverse(f), y)
        @test structisapprox(x, x2)
        @test structisapprox(y, y2)
        @test isapprox(ladj_fwd, -ladj_inv, atol = 1e-10)
        @test ladj_fwd ≈ logdensityof(MeasureBase._unweighted(μ), x) - logdensityof(MeasureBase._unweighted(ν), y)

        X = rand(μ^3)
        ys = map(f, X)
        @test batchisapprox(transports(f, X), ys)
        @test batchisapprox(f.(X), ys)
        Y, ladjs = transports_with_ladj(f, X)
        @test batchisapprox(Y, ys)
        @test batchisapprox(ladjs, map(x -> with_logabsdet_jacobian(f, x)[2], X))
        @test batchisapprox(transports(inverse(f), Y), X)
        @test isempty(transports(f, X[1:0]))
    end
end

function test_smf(μ, n = 100)
    @testset "smf($μ)" begin
        # Get `n` sorted uniforms in O(n) time
        p = rand(n)
        p .+= 0:(n-1)
        p .*= inv(n)

        F(x) = smf(μ, x)
        Finv(p) = invsmf(μ, p)

        @assert issorted(p)
        x = invsmf.(μ, p)
        @test issorted(x)
        # insupportof may return a non-Bool "don't know" (NoFastInsupport):
        @test all(x_i -> insupportof(μ, x_i) != false, x)

        @test all((Finv ∘ F).(x) .≈ x)

        # The complementary and log Stieltjes measure functions, where defined:
        if csmf(μ, first(x)) isa Real
            G(x) = csmf(μ, x)
            F_inf = F(oftype(first(x), Inf))
            @test all(G.(x) .≈ F_inf .- F.(x))
            @test all(invcsmf.(μ, G.(x)) .≈ x)
            @test all(logsmf.(μ, x) .≈ log.(F.(x)))
            @test all(logcsmf.(μ, x) .≈ log.(G.(x)))
            @test all(invlogsmf.(μ, logsmf.(μ, x)) .≈ x)
            @test all(invlogcsmf.(μ, logcsmf.(μ, x)) .≈ x)
        end

        for j in 1:n
            a = rand()
            b = rand()
            a, b = minmax(a, b)
            x = Finv(a)
            y = Finv(b)
            @test μ(Interval{:open,:closed}(x, y)) ≈ (F(y) - F(x))
            # The interval forms, to the precision of the difference:
            @test isapprox(smf(μ, x, y), F(y) - F(x); atol = 64 * eps(float(max(abs(F(x)), abs(F(y))))))
            @test isapprox(smf(μ, y, x), F(x) - F(y); atol = 64 * eps(float(max(abs(F(x)), abs(F(y))))))
            @test invsmf(μ, x, smf(μ, x, y)) ≈ y && invsmf(μ, y, smf(μ, y, x)) ≈ x
        end
    end
end

end # module MeasureBaseTestExt
