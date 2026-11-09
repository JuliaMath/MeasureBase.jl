# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test
using MeasureBase
using MeasureBase: GenContext, StdMeasure
using MeasureBase.Interface: test_interface, test_transport
using ChangesOfVariables: with_logabsdet_jacobian
using InverseFunctions: inverse
using LinearAlgebra: normalize
using LogExpFunctions: log1pexp
using Random: Xoshiro
using Statistics: std

# Measure types with exactly the methods that the extension recipes on the
# concepts page list.

# A scalar-variate leaf, the Laplace distribution, with its standard
# transport and a tail-accurate direct transport to StdLogistic:
struct _RecipeLaplace <: AbstractMeasure end
MeasureBase.logdensity_def(::_RecipeLaplace, x) = -abs(x)
MeasureBase.basemeasure(::_RecipeLaplace) = weightedmeasure(-log(2.0), Lebesgue())
MeasureBase.insupportof(::_RecipeLaplace, x) = true
MeasureBase.massof(::_RecipeLaplace) = MeasureBase.static(1.0)
MeasureBase.mspace_flatsize(::Type{<:_RecipeLaplace}) = ()
MeasureBase.preferred_stdmeasure(::Type{<:_RecipeLaplace}) = StdUniform
MeasureBase.transport_to_std(::Type{StdUniform}, ::_RecipeLaplace, x) = x < 0 ? exp(x) / 2 : 1 - exp(-x) / 2
MeasureBase.transport_from_std(::Type{StdUniform}, ::_RecipeLaplace, u) = u < 0.5 ? log(2u) : -log(2 - 2u)

_laplace_to_logistic(x) = x < 0 ? x - log(2 - exp(x)) : x + log(2 - exp(-x))
_logistic_to_laplace(y) = y < 0 ? log(2.0) - log1pexp(-y) : log1pexp(y) - log(2.0)
MeasureBase.transport_def(::StdLogistic, ::_RecipeLaplace, x) = _laplace_to_logistic(x)
MeasureBase.transport_def(::_RecipeLaplace, ::StdLogistic, y) = _logistic_to_laplace(y)
MeasureBase.batched_transport_def(ν::StdLogistic, μ::_RecipeLaplace, X) = broadcast(x -> MeasureBase.transport_def(ν, μ, x), X)
MeasureBase.batched_transport_def(ν::_RecipeLaplace, μ::StdLogistic, Y) = broadcast(y -> MeasureBase.transport_def(ν, μ, y), Y)

# An array-variate leaf, a normal distribution with independent components
# of scales σ, with its own random variate generation. The batched hooks
# keep the argument annotations of the generic methods (`Z::AbstractArray`,
# `sz::MeasureBase.SizeLike`) and so aren't ambiguous with them:
struct _RecipeDiagNormal{V<:AbstractVector{<:Real}} <: AbstractMeasure
    σ::V
end
_diagnormal_logc(d) = -sum(log, d.σ) - length(d.σ) * log(2π) / 2
_diagnormal_ld(d, X) = sum(abs2, X ./ d.σ; dims = 1) ./ -2 .+ _diagnormal_logc(d)
MeasureBase.mspace_ndims(::Type{<:_RecipeDiagNormal}) = 1
MeasureBase.mspace_flatsize(d::_RecipeDiagNormal) = size(d.σ)
MeasureBase.basemeasure(d::_RecipeDiagNormal) = weightedmeasure(_diagnormal_logc(d), Lebesgue()^length(d.σ))
MeasureBase.insupportof(::_RecipeDiagNormal, x) = true
MeasureBase.massof(::_RecipeDiagNormal) = MeasureBase.static(1.0)
MeasureBase.logdensity_def(d::_RecipeDiagNormal, x) = -sum(abs2, x ./ d.σ) / 2
MeasureBase.batched_logdensityof_impl(d::_RecipeDiagNormal, x::AbstractVector) = only(_diagnormal_ld(d, x))
MeasureBase.batched_logdensityof_impl(d::_RecipeDiagNormal, X::AbstractArray) = dropdims(_diagnormal_ld(d, X); dims = 1)
MeasureBase.preferred_stdmeasure(::Type{<:_RecipeDiagNormal}) = StdNormal
MeasureBase.transport_to_std(::Type{StdNormal}, d::_RecipeDiagNormal, x) = x ./ d.σ
MeasureBase.transport_from_std(::Type{StdNormal}, d::_RecipeDiagNormal, z) = d.σ .* z
MeasureBase.batched_transport_to_std(::Type{StdNormal}, d::_RecipeDiagNormal, X) = X ./ d.σ
MeasureBase.batched_transport_from_std(::Type{StdNormal}, d::_RecipeDiagNormal, Z::AbstractArray) = d.σ .* Z
function MeasureBase.batched_rand_impl(ctx::GenContext, d::_RecipeDiagNormal, sz::MeasureBase.SizeLike)
    d.σ .* MeasureBase.batched_rand_impl(ctx, StdNormal()^length(d.σ), sz)
end

# A measure with variates of value-dependent size, x = [a; b] with
# a ~ StdExponential() and b ~ StdNormal()^k for k = max(1, ceil(Int, a)):
struct _RecipeChain <: AbstractMeasure end
_chain_tail(a) = StdNormal()^max(1, ceil(Int, a))
_consumed(x, x_rest) = view(x, firstindex(x):(lastindex(x) - length(x_rest)))

function MeasureBase.logdensityof_with_rest(::_RecipeChain, x::AbstractVector)
    ℓ_a, a, x2 = MeasureBase.logdensityof_with_rest(StdExponential(), x)
    ℓ_b, _, x_rest = MeasureBase.logdensityof_with_rest(_chain_tail(only(a)), x2)
    return ℓ_a + ℓ_b, _consumed(x, x_rest), x_rest
end
function MeasureBase.transport_to_std_with_rest(::Type{S}, ::_RecipeChain, x::AbstractVector) where {S<:StdMeasure}
    z_a, a, x2 = MeasureBase.transport_to_std_with_rest(S, StdExponential(), x)
    z_b, _, x_rest = MeasureBase.transport_to_std_with_rest(S, _chain_tail(only(a)), x2)
    return vcat(z_a, z_b), _consumed(x, x_rest), x_rest
end
function MeasureBase.transport_from_std_with_rest(::Type{S}, ::_RecipeChain, z::AbstractVector) where {S<:StdMeasure}
    a, z2 = MeasureBase.transport_from_std_with_rest(S, StdExponential(), z)
    b, z_rest = MeasureBase.transport_from_std_with_rest(S, _chain_tail(only(a)), z2)
    return vcat(a, b), z_rest
end
function MeasureBase.batched_logdensityof_with_rest(μ::_RecipeChain, x::AbstractVector, ::Tuple{})
    ℓ, _, x_rest = MeasureBase.logdensityof_with_rest(μ, x)
    return ℓ, x_rest
end
function MeasureBase.batched_transport_to_std_with_rest(::Type{S}, μ::_RecipeChain, x::AbstractVector, ::Tuple{}) where {S<:StdMeasure}
    z, _, x_rest = MeasureBase.transport_to_std_with_rest(S, μ, x)
    return z, x_rest
end
function MeasureBase.batched_transport_from_std_with_rest(::Type{S}, μ::_RecipeChain, z::AbstractVector, ::Tuple{}) where {S<:StdMeasure}
    MeasureBase.transport_from_std_with_rest(S, μ, z)
end
MeasureBase.fixed_stream_size(::Type{<:_RecipeChain}) = MeasureBase.static(false)

MeasureBase.insupportof(::_RecipeChain, x) = MeasureBase.NoFastInsupport{_RecipeChain}()
MeasureBase.massof(::_RecipeChain) = MeasureBase.static(1.0)
MeasureBase.preferred_stdmeasure(::Type{<:_RecipeChain}) = StdNormal
function MeasureBase.rand_impl(ctx::GenContext, ::_RecipeChain)
    a = MeasureBase.rand_impl(ctx, StdExponential())
    return vcat(a, MeasureBase.rand_impl(ctx, _chain_tail(a)))
end
MeasureBase.mspace_ndims(::Type{<:_RecipeChain}) = 1


@testset "extension recipes" begin
    @testset "scalar-variate leaf" begin
        L = _RecipeLaplace()
        ℓ(x) = -abs(x) - log(2.0)
        xs = [-1.5, 0.3, 2.0]
        @test logdensityof(L, 0.3) ≈ ℓ(0.3)
        @test logdensities(L, xs) ≈ ℓ.(xs)
        @test logdensityof(L^3, xs) ≈ sum(ℓ, xs)
        @test logdensityof(productmeasure((L, StdNormal())), (0.3, 0.1)) ≈ ℓ(0.3) + logdensityof(StdNormal(), 0.1)
        @test logdensityof(mcombine(vcat, L, StdNormal()^2), [0.3, 0.1, 0.2]) ≈ ℓ(0.3) + logdensityof(StdNormal()^2, [0.1, 0.2])
        @test normalize(L) === L
        test_interface(L)
        test_transport(StdNormal(), L)
        test_transport(L, StdNormal())
        # The pushforward of a measure by its transport to another measure
        # is that measure:
        f = transport_to(StdNormal(), L)
        @test logdensityof(pushfwd(f, L), 0.7) ≈ logdensityof(StdNormal(), 0.7)
    end

    @testset "direct transports between pairs of measure types" begin
        L = _RecipeLaplace()
        f = transport_to(StdLogistic(), L)
        X = [-0.5, 40.0, 100.0]
        # The direct transport stays tail-accurate, the one via the
        # StdUniform pivot saturates:
        @test f(100.0) ≈ 100.0 + log(2 - exp(-100.0))
        @test transports(f, X) == f.(X) == map(f, X)
        Y, ladjs = transports_with_ladj(f, X)
        @test Y == map(f, X)
        @test ladjs ≈ log(2.0) .- log.(2 .- exp.(-abs.(X)))
        @test ladjs ≈ map(x -> with_logabsdet_jacobian(f, x)[2], X)
        @test transports(inverse(f), Y) ≈ X
        test_transport(StdLogistic(), L)
        test_transport(L, StdLogistic())
    end

    @testset "array-variate leaf" begin
        σ = [0.5, 1.0, 2.0]
        d = _RecipeDiagNormal(σ)
        x, y = [0.1, 0.2, -0.3], [1.0, -0.5, 0.2]
        ℓ(x) = sum(logdensityof.(StdNormal(), x ./ σ) .- log.(σ))
        @test logdensityof(d, x) ≈ ℓ(x)
        @test logdensities(d, hcat(x, y)) ≈ [ℓ(x), ℓ(y)]
        @test logdensityof(d^2, hcat(x, y)) ≈ ℓ(x) + ℓ(y)
        @test logdensity_rel(d, StdNormal()^3, x) ≈ ℓ(x) - logdensityof(StdNormal()^3, x)
        @test logdensityof(productmeasure((d, StdNormal())), (x, 0.5)) ≈ ℓ(x) + logdensityof(StdNormal(), 0.5)
        @test logdensityof(mcombine(vcat, d, StdNormal()), [x; 0.5]) ≈ ℓ(x) + logdensityof(StdNormal(), 0.5)
        @test getdof(d) == 3 && fast_dof(d) == 3
        @test size(rand(d)) == (3,) && size(rand(d^4)) == (4,) && size(rand(d^(2, 3))) == (2, 3)
        @test size(rand(d^MeasureBase.static(2))) == (2,)
        X = MeasureBase.batched_rand_impl(GenContext{Float64}(Xoshiro(42)), d, (10^4,))
        @test size(X) == (3, 10^4)
        @test vec(std(X; dims = 2)) ≈ σ rtol = 0.05
        test_interface(d)
        test_transport(StdNormal()^3, d)
        test_transport(d, StdNormal()^3)
        test_transport(StdUniform()^3, d)
    end

    @testset "measures of value-dependent size" begin
        μ = _RecipeChain()
        x, y = [1.5, 0.3, -0.2], [0.5, 0.1]
        ℓ(x) = logdensityof(StdExponential(), x[1]) + logdensityof(StdNormal()^(length(x) - 1), x[2:end])
        @test logdensityof(μ, x) ≈ ℓ(x) && logdensityof(μ, y) ≈ ℓ(y)
        @test logdensityof(μ, [-1.0, 0.3]) == -Inf
        @test_throws ArgumentError logdensityof(μ, [0.5, 0.1, 0.2])
        @test logdensities(μ, [x, y]) ≈ [ℓ(x), ℓ(y)]
        @test logdensityof(μ^2, [x, y]) ≈ ℓ(x) + ℓ(y)
        ν = mcombine(vcat, μ, StdNormal())
        ℓ_ν(x) = ℓ(x[1:(end - 1)]) + logdensityof(StdNormal(), x[end])
        @test logdensities(ν, [[x; 0.7], [y; 0.7]]) ≈ [ℓ_ν([x; 0.7]), ℓ_ν([y; 0.7])]
        @test logdensityof(mbind(_ -> StdNormal(), μ, vcat), [x; 0.7]) ≈ ℓ_ν([x; 0.7])

        f = transport_to(StdNormal()^3, μ)
        @test f(x) ≈ [transport_to(StdNormal(), StdExponential())(1.5), 0.3, -0.2]
        @test inverse(f)(f(x)) ≈ x
        @test with_logabsdet_jacobian(f, x)[2] ≈ ℓ(x) - logdensityof(StdNormal()^3, f(x))
        x2 = [1.2, -0.4, 0.6]
        @test transports(f, [x, x2]) ≈ [f(x), f(x2)]
        @test transport_to(StdNormal()^4, ν)([x; 0.7]) ≈ [f(x); 0.7]

        # Without a base measure, densities, transports, the degrees of
        # freedom and the argument checks come from the with-rest hooks:
        @test getdof(μ) isa MeasureBase.NoDOF && fast_dof(μ) isa MeasureBase.NoDOF
        xs = [rand(μ) for _ in 1:5]
        @test all(x -> length(x) == 1 + max(1, ceil(Int, x[1])), xs)
        @test all(isfinite, logdensities(μ, xs))
    end
end

