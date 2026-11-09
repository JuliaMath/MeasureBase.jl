using Test
using MeasureBase
using MeasureBase: LebesgueBase
using InverseFunctions: inverse
using LinearAlgebra: I
using Statistics: mean, var
import AffineMaps
import Distributions

# A scalar-variate leaf with exactly the declarations the concepts page asks
# for, the Laplace distribution with its standard transport:
struct _MinimalLaplace <: AbstractMeasure end
MeasureBase.logdensity_def(::_MinimalLaplace, x) = -abs(x)
MeasureBase.basemeasure(::_MinimalLaplace) = weightedmeasure(-log(2.0), LebesgueBase())
MeasureBase.insupportof(::_MinimalLaplace, x) = true
MeasureBase.mspace_flatsize(::Type{<:_MinimalLaplace}) = ()
MeasureBase.massof(::_MinimalLaplace) = MeasureBase.static(1.0)
MeasureBase.preferred_stdmeasure(::Type{<:_MinimalLaplace}) = StdUniform
MeasureBase.transport_to_std(::Type{StdUniform}, ::_MinimalLaplace, x) = x < 0 ? exp(x) / 2 : 1 - exp(-x) / 2
MeasureBase.transport_from_std(::Type{StdUniform}, ::_MinimalLaplace, u) = u < 0.5 ? log(2u) : -log(2 - 2u)
_laplace_ld(x) = -abs(x) - log(2.0)

# A Laplace leaf with a direct, tail-accurate transport to StdExponential
# that only specializes the point form `transport_def`:
struct _DirectLaplace <: AbstractMeasure end
MeasureBase.logdensity_def(::_DirectLaplace, x) = -abs(x)
MeasureBase.basemeasure(::_DirectLaplace) = weightedmeasure(-log(2.0), LebesgueBase())
MeasureBase.insupportof(::_DirectLaplace, x) = true
MeasureBase.mspace_flatsize(::Type{<:_DirectLaplace}) = ()
MeasureBase.massof(::_DirectLaplace) = MeasureBase.static(1.0)
MeasureBase.preferred_stdmeasure(::Type{<:_DirectLaplace}) = StdUniform
MeasureBase.transport_to_std(::Type{StdUniform}, ::_DirectLaplace, x) = x < 0 ? exp(x) / 2 : 1 - exp(-x) / 2
MeasureBase.transport_from_std(::Type{StdUniform}, ::_DirectLaplace, u) = u < 0.5 ? log(2u) : -log(2 - 2u)
MeasureBase.transport_def(::StdExponential, ::_DirectLaplace, x) = x > 0 ? x + log(2.0) : -log1p(-exp(x) / 2)

# The Laplace leaf without a declared variate rank:
struct _RanklessLaplace <: AbstractMeasure end
MeasureBase.logdensity_def(::_RanklessLaplace, x) = -abs(x)
MeasureBase.basemeasure(::_RanklessLaplace) = weightedmeasure(-log(2.0), LebesgueBase())
MeasureBase.insupportof(::_RanklessLaplace, x) = true
MeasureBase.preferred_stdmeasure(::Type{<:_RanklessLaplace}) = StdUniform
MeasureBase.transport_to_std(::Type{StdUniform}, ::_RanklessLaplace, x) = x < 0 ? exp(x) / 2 : 1 - exp(-x) / 2
MeasureBase.transport_from_std(::Type{StdUniform}, ::_RanklessLaplace, u) = u < 0.5 ? log(2u) : -log(2 - 2u)

# An isotropic normal leaf with vector variates that implements only the
# batched base density:
struct _BatchedIsoNormal <: AbstractMeasure end
MeasureBase.mspace_flatsize(::Type{<:_BatchedIsoNormal}) = (MeasureBase.static(2),)
MeasureBase.basemeasure(::_BatchedIsoNormal) = Lebesgue()^2
MeasureBase.insupportof(::_BatchedIsoNormal, x) = true
function MeasureBase.batched_logdensity_def(::_BatchedIsoNormal, X)
    r = sum(abs2, X; dims = 1) ./ -2 .- log(2π)
    X isa AbstractVector ? only(r) : dropdims(r; dims = 1)
end
_iso_ld(x) = -sum(abs2, x) / 2 - log(2π)

# A measure of value-dependent variate size that implements only the
# with-rest hooks: `x = [a; b]` with `a ~ StdExponential()` and
# `b ~ StdNormal()^k(a)`:
struct _ExpChain <: AbstractMeasure end
_chain_k(a) = max(1, ceil(Int, a))
MeasureBase.fixed_stream_size(::Type{<:_ExpChain}) = MeasureBase.static(false)
MeasureBase.preferred_stdmeasure(::Type{<:_ExpChain}) = StdNormal
function MeasureBase.logdensityof_with_rest(::_ExpChain, x::AbstractVector)
    ℓa, a, x2 = MeasureBase.logdensityof_with_rest(StdExponential(), x)
    ℓb, _, x_rest = MeasureBase.logdensityof_with_rest(StdNormal()^_chain_k(a), x2)
    return ℓa + ℓb, view(x, firstindex(x):(firstindex(x) + _chain_k(a))), x_rest
end
function MeasureBase.transport_to_std_with_rest(::Type{StdNormal}, ::_ExpChain, x::AbstractVector)
    za, a, x2 = MeasureBase.transport_to_std_with_rest(StdNormal, StdExponential(), x)
    zb, _, x_rest = MeasureBase.transport_to_std_with_rest(StdNormal, StdNormal()^_chain_k(a), x2)
    return vcat(za, zb), view(x, firstindex(x):(firstindex(x) + _chain_k(a))), x_rest
end
function MeasureBase.transport_from_std_with_rest(::Type{StdNormal}, ::_ExpChain, z::AbstractVector)
    a, z2 = MeasureBase.transport_from_std_with_rest(StdNormal, StdExponential(), z)
    b, z_rest = MeasureBase.transport_from_std_with_rest(StdNormal, StdNormal()^_chain_k(a), z2)
    return vcat(a, b), z_rest
end
_chain_ld(x) = -x[1] + sum(xi -> -xi^2 / 2 - log(2π) / 2, x[2:end])

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
    # Supports given as types contain only values of the type:
    for (c, out) in ((Counting(0:10), 11), (Counting(Set([1, 2, 7])), 3), (Counting([1, 2, 7]), 3), (Counting(Int), nothing))
        test_measure_contracts(c, [1, 2, 7], out; transport = false)
        test_measure_contracts(c^3, [[1, 2, 7], [7, 1, 1]], isnothing(out) ? nothing : [1, 2, out]; transport = false)
        test_measure_contracts(productmeasure((c, StdNormal())), [(1, 0.5), (7, -0.5)], isnothing(out) ? nothing : (out, 0.5); transport = false)
        test_composition_invariants(c, 1, 2)
    end
    @test logdensityof(Counting(Int), 1.5) == logdensityof(Counting(Int)^2, [1, 1.5]) == -Inf
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
        @test MeasureBase.with_logabsdet_jacobian(f, 0.3)[2] ≈ _laplace_ld(0.3) - logdensityof(StdNormal(), f(0.3))
    end

    @testset "array leaves with only a batched base density" begin
        d = _BatchedIsoNormal()
        x, y = [0.1, 0.2], [-1.0, 0.5]
        @test logdensityof(d, x) ≈ MeasureBase.logdensity_def(d, x) ≈ _iso_ld(x)
        @test logdensities(d, hcat(x, y)) ≈ [_iso_ld(x), _iso_ld(y)]
        @test logdensityof(d^2, hcat(x, y)) ≈ _iso_ld(x) + _iso_ld(y)
        @test logdensity_rel(d, Lebesgue()^2, x) ≈ _iso_ld(x)
        @test logdensityof(superpose(d, d), x) ≈ _iso_ld(x) + log(2)
        test_composition_invariants(d, x, y)
    end

    @testset "measures that implement only the stream protocol" begin
        μ = _ExpChain()
        x1, x2 = [0.5, 0.3], [1.5, 0.3, -0.2]
        @test logdensityof(μ, x1) ≈ _chain_ld(x1)
        @test logdensityof(μ, x2) ≈ _chain_ld(x2)
        @test_throws ArgumentError logdensityof(μ, [0.5, 0.3, 0.1])
        @test logdensities(μ, [x1, x2]) ≈ _chain_ld.([x1, x2])
        @test logdensities(μ, hcat(x1, [0.7, -0.3])) ≈ [_chain_ld(x1), _chain_ld([0.7, -0.3])]
        @test logdensityof(mcombine(vcat, μ, StdNormal()), vcat(x2, 0.4)) ≈ _chain_ld(x2) + logdensityof(StdNormal(), 0.4)
        @test logdensityof(MeasureBase.mbind(_ -> μ, StdNormal(), vcat), vcat(0.1, x2)) ≈ logdensityof(StdNormal(), 0.1) + _chain_ld(x2)
        @test logdensityof(μ^2, [x1, x2]) ≈ _chain_ld(x1) + _chain_ld(x2)
        @test MeasureBase.getdof(μ) isa MeasureBase.NoDOF && MeasureBase.fast_dof(μ) isa MeasureBase.NoDOF
        # The first standard variate is the normal quantile of the exponential cdf:
        f = transport_to(StdNormal()^3, μ)
        z = f(x2)
        @test z[1] ≈ Distributions.quantile(Distributions.Normal(), 1 - exp(-1.5)) && z[2:3] == x2[2:3]
        @test inverse(f)(z) ≈ x2
    end

    @testset "leaves without a declared variate rank" begin
        L = _RanklessLaplace()
        @test logdensityof(L, 0.3) ≈ _laplace_ld(0.3)
        @test transport_to(StdUniform(), L)(0.3) ≈ 1 - exp(-0.3) / 2
        # Draws and transports from standard measures need the rank:
        @test_throws ArgumentError rand(L)
        @test_throws ArgumentError testvalue(L)
        @test_throws ArgumentError rand(L^2)
        @test_throws ArgumentError transport_to(L, StdUniform())(0.3)
    end

    @testset "batches follow direct point transports" begin
        f = transport_to(StdExponential(), _DirectLaplace())
        X = [-0.7, 0.5, 40.0, 100.0]
        # The exponential quantile of the Laplace cdf, exact in the upper tail:
        @test map(f, X) ≈ [-log1p(-exp(-0.7) / 2), 0.5 + log(2), 40 + log(2), 100 + log(2)]
        @test transports(f, X) == f.(X) == map(f, X)
        Y, ladjs = transports_with_ladj(f, X)
        @test Y == map(f, X)
        @test ladjs ≈ map(x -> MeasureBase.with_logabsdet_jacobian(f, x)[2], X)
        @test ladjs ≈ logdensityof.(_DirectLaplace(), X) .- logdensityof.(StdExponential(), Y)
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
