# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

export logdensities

"""
    logdensities(μ::AbstractMeasure, X::AbstractArray)

Compute the log-density of `μ` at each point in `X`.

Returns an array of the shape of `X`, semantically equivalent to
`logdensityof.(Ref(μ), X)`. The computation may be fused across points,
though: power measures with flat variate storage (e.g. based on
`ArraysOfArrays.ArrayOfSimilarArrays`) evaluate as a single flat broadcast
plus a segmented reduction over the underlying flat data, compatible with
GPU-backed storage.

Measure types should specialize [`MeasureBase.batched_logdensityof_impl`](@ref)
instead of `logdensities` itself.
"""
function logdensities end

@inline function logdensities(μ::AbstractMeasure, X::AbstractArray)
    _logdensities(logdensityof_impl, μ, X)
end

"""
    MeasureBase.batched_logdensityof_impl(μ::AbstractMeasure, X::AbstractArray)

Implements [`logdensities(μ, X)`](@ref logdensities) for arrays `X` of
plain `μ`-variates. Power measures never reach `batched_logdensityof_impl`, their
power structure is processed generically beforehand.

Measure types that support fused multi-point evaluation should specialize
`batched_logdensityof_impl`. Implementations must preserve the shape of `X` and
must handle points outside the support of `μ` (the result must be `-Inf`
at such points).

The default implementation broadcasts the log-density over `X` for
measures with scalar variates and falls back to a `map` over `X`
otherwise.
"""
function batched_logdensityof_impl end

function batched_logdensityof_impl(μ::AbstractMeasure, X::AbstractArray)
    _logdensities_generic(logdensityof_impl, μ, X)
end

# Batched density machinery, parameterized over the point-level density
# function `f` (`logdensityof_impl` or `logdensity_def`).
#
# `_logdensities(f, μ, X, powers...)` treats each element of `X` as a
# variate of `(μ^pN)^…^p1` for `powers = (p1, …, pN)` (power axes ordered
# outermost first, i.e. in the order in which they are encountered when
# descending into a variate) and returns the density sum for each element,
# preserving the shape of `X`. Power measures are unwrapped into the power
# axes arguments before any other dispatch happens, so implementation
# methods only ever dispatch on plain measure types.

@inline function _logdensities(f::F, μ, X::AbstractArray, powers::Vararg{Any,N}) where {F,N}
    _logdensities_stripped(f, μ, X, powers...)
end

@inline function _logdensities(
    f::F,
    μ::PowerMeasure,
    X::AbstractArray,
    powers::Vararg{Any,N},
) where {F,N}
    _logdensities(f, pwr_base(μ), X, powers..., pwr_axes(μ))
end

@inline function _logdensities_stripped(f::F, μ, X::AbstractArray) where {F}
    _batched_logdensityof_impl(f, μ, X)
end

function _logdensities_stripped(
    f::F,
    μ,
    X::AbstractArray,
    p1,
    powers::Vararg{Any,N},
) where {F,N}
    map(x -> _powered_ld(f, μ, x, p1, powers...), X)
end

function _logdensities_stripped(
    f::F,
    μ,
    X::ArrayOfSimilarArrays{<:Real},
    p1,
    powers::Vararg{Any,N},
) where {F,N}
    _logdensities_fused(f, μ, X, mspace_elsize(μ), p1, powers...)
end

# Absolute densities go through the `batched_logdensityof_impl` extension point:
@inline function _batched_logdensityof_impl(::typeof(logdensityof_impl), μ, X::AbstractArray)
    batched_logdensityof_impl(μ, X)
end

@inline function _batched_logdensityof_impl(f::F, μ, X::AbstractArray) where {F}
    _logdensities_generic(f, μ, X)
end

@inline function _logdensities_generic(f::F, μ, X::AbstractArray) where {F}
    _logdensities_byelsize(f, μ, X, mspace_elsize(μ))
end

# Scalar variates evaluate as a single flat broadcast:
@inline function _logdensities_byelsize(
    f::F,
    μ,
    X::AbstractArray{<:Real},
    ::Tuple{},
) where {F}
    broadcast(Base.Fix1(f, μ), X)
end

@inline function _logdensities_byelsize(f::F, μ, X::AbstractArray, ::Any) where {F}
    map(Base.Fix1(f, μ), X)
end

# Scalar-variate measure with flat variate storage: evaluate as a single
# flat broadcast followed by a segmented reduction over the per-point
# power structure:
function _logdensities_fused(
    f::F,
    μ,
    X::ArrayOfSimilarArrays{<:Real,M},
    ::Tuple{},
    powers::Vararg{Any,N},
) where {F,M,N}
    sz_inner = _flat_powers_size(powers...)
    if length(sz_inner) == M
        X_flat = flatview(X)
        if ntuple(i -> size(X_flat, i), Val(M)) != sz_inner
            throw(ArgumentError("Size of variates doesn't match size of power measure"))
        end
        ld_flat = broadcast(Base.Fix1(f, μ), X_flat)
        reshape(sum(ld_flat, dims = ntuple(identity, Val(M))), size(X))
    else
        map(x -> _powered_ld(f, μ, x, powers...), X)
    end
end

function _logdensities_fused(
    f::F,
    μ,
    X::AbstractArray,
    ::Any,
    p1,
    powers::Vararg{Any,N},
) where {F,N}
    map(x -> _powered_ld(f, μ, x, p1, powers...), X)
end

# The flat size of a variate of `(μ^pN)^…^p1` for a scalar-variate `μ`,
# innermost power axes vary fastest:
@inline _flat_powers_size() = ()
@inline function _flat_powers_size(p1, powers::Vararg{Any,N}) where {N}
    (_flat_powers_size(powers...)..., axes2size(p1)...)
end

# Scalar counterpart of `_logdensities`: log-density of `(μ^pN)^…^p1` at a
# single variate `x`.
@inline _powered_ld(f::F, μ, x) where {F} = f(μ, x)

@inline function _powered_ld(
    f::F,
    μ::PowerMeasure,
    x,
    p1,
    powers::Vararg{Any,N},
) where {F,N}
    _powered_ld(f, pwr_base(μ), x, p1, powers..., pwr_axes(μ))
end

@inline function _powered_ld(f::F, μ, x, p1, powers::Vararg{Any,N}) where {F,N}
    if axes2size(p1) != maybestatic_size(x)
        throw(ArgumentError("Size of variate doesn't match size of power measure"))
    end
    R = _powered_ld_type(f, μ, x, powers...)
    if isempty(x)
        zero(R)::R
    else
        # Conversion needed since summation can turn static into dynamic values:
        convert(R, _powered_ld_sum(f, μ, x, powers...))::R
    end
end

@inline _powered_ld_sum(f::F, μ, x) where {F} = sum(Base.Fix1(f, μ), x)

@inline function _powered_ld_sum(f::F, μ, x, p1, powers::Vararg{Any,N}) where {F,N}
    sum(_logdensities(f, μ, x, p1, powers...))
end

@inline function _powered_ld_type(f::F, ::MU, x, powers::Vararg{Any,N}) where {F,MU,N}
    Core.Compiler.return_type(_powered_ld, Tuple{F,MU,eltype(x),map(typeof, powers)...})
end
