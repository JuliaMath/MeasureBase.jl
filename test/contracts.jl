using Test
using MeasureBase
using MeasureBase: LebesgueBase
using InverseFunctions: inverse
using LinearAlgebra: I
using Statistics: mean, var
import AffineMaps

# A scalar-variate leaf with exactly the declarations the concepts page asks
# for, the Laplace distribution with its standard transport:
struct _MinimalLaplace <: AbstractMeasure end
MeasureBase.logdensity_def(::_MinimalLaplace, x) = -abs(x)
MeasureBase.basemeasure(::_MinimalLaplace) = weightedmeasure(-log(2.0), LebesgueBase())
MeasureBase.insupportof(::_MinimalLaplace, x) = true
MeasureBase.mspace_flatsize(::Type{<:_MinimalLaplace}) = ()
MeasureBase.mspace_flatsize(::_MinimalLaplace) = ()
MeasureBase.preferred_stdmeasure(::Type{<:_MinimalLaplace}) = StdUniform
MeasureBase.transport_to_std(::Type{StdUniform}, ::_MinimalLaplace, x) = x < 0 ? exp(x) / 2 : 1 - exp(-x) / 2
MeasureBase.transport_from_std(::Type{StdUniform}, ::_MinimalLaplace, u) = u < 0.5 ? log(2u) : -log(2 - 2u)
_laplace_ld(x) = -abs(x) - log(2.0)

# Variates of `μ^2` from two variates of `μ`, and batches of two of those:
_pwr2(x1::Real, x2::Real) = [x1, x2]
_pwr2(x1::AbstractVector, x2::AbstractVector) = hcat(x1, x2)
_batch2(X, Y) = cat(X, Y; dims = ndims(X) + 1)

# Densities through powers and base measures: densities of powers are sums
# over their elements, densities relative to the base measure and the
# density of the base add up to the density, and so do relative densities
# of powers relative to the powers of the base measures.
function test_composition_invariants(μ, x1, x2)
    @testset "composition invariants of $μ" begin
        ν = basemeasure(μ)
        ℓ1, ℓ2 = logdensityof(μ, x1), logdensityof(μ, x2)
        d1, d2 = logdensity_def(μ, x1), logdensity_def(μ, x2)
        X, Y = _pwr2(x1, x2), _pwr2(x2, x1)
        @test ℓ1 ≈ d1 + logdensityof(ν, x1)
        @test logdensityof(μ^2, X) ≈ ℓ1 + ℓ2
        @test logdensity_def(μ^2, X) ≈ d1 + d2
        @test logdensity_rel(μ^2, ν^2, X) ≈ d1 + d2
        @test logdensities(μ^2, _batch2(X, Y)) ≈ [ℓ1 + ℓ2, ℓ1 + ℓ2]
        @test MeasureBase.batched_logdensity_def(μ^2, _batch2(X, Y)) ≈ [d1 + d2, d1 + d2]
    end
end

# Checks the batching, support and precision contracts of measure `μ` on
# the variates `xs`, `x_out` (if given) must lie outside of the support:
function test_measure_contracts(μ, xs::AbstractVector, x_out; transport::Bool = true)
    @testset "contracts of $μ" begin
        _f32(x::Real) = Float32(x)
        _f32(x::AbstractArray) = Float32.(x)
        _f32(x::Union{Tuple,NamedTuple}) = map(_f32, x)
        xs32 = map(_f32, xs)

        ld = map(x -> logdensityof(μ, x), xs)
        @test all(isfinite, ld)
        @test logdensities(μ, xs) ≈ ld
        @test isempty(logdensities(μ, xs[1:0]))
        @test logdensityof(μ, first(xs32)) isa Float32
        @test eltype(logdensities(μ, xs32)) == Float32
        if !isnothing(x_out)
            @test logdensityof(μ, x_out) == -Inf
            @test logdensities(μ, [x_out, first(xs)]) ≈ [-Inf, first(ld)]
        end

        if transport
            f = transport_to(StdUniform, μ)
            ys = map(f, xs)
            @test all(map((a, b) -> all(a .≈ b), f.(xs), ys))
            @test all(map((a, b) -> all(map(≈, a, b)), inverse(f).(ys), xs))
            @test isempty(f.(xs[1:0]))
            @test eltype(f(first(xs32))) == Float32
            isnothing(x_out) || @test any(isnan, f(x_out))
        end
    end
end

test_measure_contracts(μ, x_out; kwargs...) =
    test_measure_contracts(μ, [rand(μ) for _ in 1:5], x_out; kwargs...)

@testset "measure contracts" begin
    test_measure_contracts(StdNormal(), nothing)
    test_measure_contracts(StdLogistic(), nothing)
    test_measure_contracts(StdExponential(), -1.0)
    test_measure_contracts(StdUniform(), 1.5)
    test_measure_contracts(StdExponential()^3, [1.0, -1.0, 1.0])
    test_measure_contracts(StdExponential()^(2, 2), [1.0 -1.0; 1.0 1.0])
    test_measure_contracts(weightedmeasure(0.5, StdExponential()), -1.0)
    test_measure_contracts(half(StdNormal()), -1.0)
    test_measure_contracts(pushfwd(exp, StdNormal()), nothing)
    test_measure_contracts(productmeasure((StdNormal(), StdExponential())), (0.5, -1.0))
    test_measure_contracts(productmeasure((a = StdNormal(), b = StdUniform())), (a = 0.5, b = 2.0))
    test_measure_contracts(mcombine(vcat, StdNormal()^2, StdExponential()^2), [0.1, 0.2, -1.0, 1.0])
    test_measure_contracts(Dirac(2.0), 3.0; transport = false)
    test_measure_contracts(superpose(StdUniform(), StdExponential()), -1.0; transport = false)
    test_measure_contracts(restrict(>(0), StdNormal()), [0.5, 1.0, 2.0], -1.0; transport = false)
    test_measure_contracts(mintegrate_exp(x -> -x, StdExponential()), [0.5, 1.0, 2.0], -1.0; transport = false)

    test_composition_invariants(StdNormal(), 0.3, -1.2)
    test_composition_invariants(weightedmeasure(0.5, StdNormal()), 0.3, -1.2)
    test_composition_invariants(half(StdNormal()), 0.3, 1.2)
    test_composition_invariants(restrict(>(0), StdExponential()), 0.3, 1.2)
    test_composition_invariants(mintegrate_exp(x -> -x, StdExponential()), 0.3, 1.2)
    test_composition_invariants(pushfwd(exp, StdNormal()), 0.3, 1.2)
    test_composition_invariants(pushfwd(AffineMaps.Mul(2.0), StdNormal()), 0.3, -1.2)
    test_composition_invariants(pushfwd(AffineMaps.Mul(2.0 * I(2)), StdNormal()^2), [0.3, -1.2], [0.5, 0.1])
    test_composition_invariants(StdNormal()^3, [0.3, -1.2, 0.5], [0.5, 0.1, 2.0])
    test_composition_invariants(_MinimalLaplace(), 0.3, -1.2)

    @testset "relative densities of measures of equal type" begin
        a, b = restrict(>(0), StdNormal()), restrict(>(1), StdNormal())
        @test typeof(a) == typeof(b)
        @test logdensity_rel(a, b, 2.0) == 0
        @test logdensity_rel(a, b, 0.5) == Inf
        @test logdensity_rel(b, a, 0.5) == -Inf
        @test isnan(logdensity_rel(a, b, -1.0))
        @test logdensity_rel(Dirac(0.0), Dirac(1.0), 0.0) == Inf
        @test logdensity_rel(Dirac(0.0), Dirac(1.0), 1.0) == -Inf
        @test logdensity_rel(2.0 * StdNormal(), 3.0 * StdNormal(), 0.7) ≈ log(2 / 3)
        p2, p3 = pushfwd(AffineMaps.Mul(2.0), StdNormal()), pushfwd(AffineMaps.Mul(3.0), StdNormal())
        ℓN(x) = logdensityof(StdNormal(), x)
        @test logdensity_rel(p2, p3, 0.7) ≈ (ℓN(0.7 / 2) - log(2)) - (ℓN(0.7 / 3) - log(3))
    end

    @testset "measures with only the documented minimal declarations" begin
        L = _MinimalLaplace()
        xs = [-1.5, 0.3, 2.0]
        @test logdensityof(L, 0.3) ≈ _laplace_ld(0.3)
        @test logdensities(L, xs) ≈ _laplace_ld.(xs)
        @test logdensityof(productmeasure((L, StdNormal())), (0.3, 0.1)) ≈ _laplace_ld(0.3) + logdensityof(StdNormal(), 0.1)
        @test logdensities(productmeasure((L, StdNormal())), (xs, xs)) ≈ _laplace_ld.(xs) .+ logdensityof.(StdNormal(), xs)
        @test logdensityof(mcombine(vcat, L, StdNormal()^2), [0.3, 0.1, 0.2]) ≈ _laplace_ld(0.3) + logdensityof(StdNormal()^2, [0.1, 0.2])
        X = rand(L^(10^4))
        @test isapprox(mean(X), 0, atol = 0.1) && isapprox(var(X), 2, rtol = 0.1)
        f = transport_to(StdNormal(), L)
        @test inverse(f)(f(0.3)) ≈ 0.3
        @test f.(xs) ≈ map(f, xs)
    end

    # Pushforwards are only defined on the image of their map, the inverse
    # may throw outside of it:
    @test logdensityof(pushfwd(exp, StdNormal()), 2.0) ≈ -log(2.0) - log(2π) / 2 - log(2.0)^2 / 2
    @test_throws DomainError logdensityof(pushfwd(exp, StdNormal()), -1.0)

    # Known gaps:
    let P = productmeasure([StdNormal(), StdExponential()])
        @test isempty(logdensities(P, rand(P^0)))
    end
end
