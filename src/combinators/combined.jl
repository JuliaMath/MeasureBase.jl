"""
    MeasureBase.tpmeasure_split_combined(f_c, α::AbstractMeasure, ab)

Splits a combined value `ab` that originated from combining a point `a`
from the space of a measure `α` with a point `b` from the space of
another measure `β` via `ab = f_c(a, b)`.

Returns a semantic equivalent of
`(MeasureBase.transportmeasure(α, a), a, b)`.

With `a_orig = rand(α)`, `b_orig = rand(β)` and
`ab = f_c(a_orig, b_orig)`, the following must hold true:

```julia
tpm_α, a, b = tpmeasure_split_combined(f_c, α, ab)
a ≈ a_orig && b ≈ b_orig
```
"""
function tpmeasure_split_combined end

function tpmeasure_split_combined(f_c, α::AbstractMeasure, ab)
    a, b = _generic_split_combined(f_c, α, ab)
    return transportmeasure(α, a), a, b
end

@inline _generic_split_combined(::typeof(tuple), ::AbstractMeasure, x::Tuple{Vararg{Any,2}}) = x
@inline _generic_split_combined(::Type{Pair}, ::AbstractMeasure, ab::Pair) = (ab...,)

function _generic_split_combined(f_c::FC, α::AbstractMeasure, ab) where {FC}
    _split_variate_byvalue(f_c, testvalue(α), ab)
end

_split_variate_byvalue(::typeof(vcat), test_a::AbstractVector, ab::AbstractVector) =
    _split_after(ab, length(test_a))

_split_variate_byvalue(::typeof(vcat), ::Number, ab::AbstractVector) =
    _consume_from_stream(ab, ())

_split_variate_byvalue(::typeof(vcat), ::NTuple{N,Any}, ab::Tuple) where {N} =
    _split_after(ab, Val{N}())

function _split_variate_byvalue(::typeof(merge), ::NamedTuple{names_a}, ab::NamedTuple) where {names_a}
    _split_after(ab, Val(names_a))
end


@doc raw"""
    mcombine(f_c, α::AbstractMeasure, β::AbstractMeasure)

Combines two measures `α` and `β` to a combined measure via a point
combination function `f_c`.

`f_c` must combine a given point `a` from the space of measure `α` with a
given point `b` from the space of measure `β` to a single value
`ab = f_c(a, b)` in the space of the combined measure
`μ = mcombine(f_c, α, β)`.

The combined measure has the mathematical interpretation (on sets
$$A$$ and $$B$$)

```math
\mu(f_c(A, B)) = \alpha(A)\, \beta(B)
```
"""
function mcombine end
export mcombine

@inline function mcombine(f_c, α::AbstractMeasure, β::AbstractMeasure)
    _generic_mcombine_impl_stage1(f_c, α, β)
end

@inline _generic_mcombine_impl_stage1(::typeof(firstarg), α::AbstractMeasure, β::AbstractMeasure) = α
@inline _generic_mcombine_impl_stage1(::typeof(secondarg), α::AbstractMeasure, β::AbstractMeasure) = β

@inline function _generic_mcombine_impl_stage1(::typeof(tuple), α::AbstractMeasure, β::AbstractMeasure)
    productmeasure((α, β))
end

@inline function _generic_mcombine_impl_stage1(
    f_c::Union{typeof(vcat),typeof(merge)},
    α::AbstractProductMeasure,
    β::AbstractProductMeasure,
)
    _mcombine_product_shortcut(f_c, marginals(α), marginals(β), α, β)
end

function _mcombine_product_shortcut(::typeof(vcat), ma::AbstractVector{T}, mb::AbstractVector{T}, α, β) where {T}
    isconcretetype(T) ? productmeasure(vcat(ma, mb)) : _generic_mcombine_impl_stage2(vcat, α, β)
end
_mcombine_product_shortcut(::typeof(merge), ma::NamedTuple, mb::NamedTuple, α, β) =
    productmeasure(merge(ma, mb))
_mcombine_product_shortcut(f_c, ma, mb, α, β) = _generic_mcombine_impl_stage2(f_c, α, β)

@inline function _generic_mcombine_impl_stage1(f_c, α::AbstractMeasure, β::AbstractMeasure)
    _generic_mcombine_impl_stage2(f_c, α, β)
end

@inline function _generic_mcombine_impl_stage2(f_c, α::AbstractMeasure, β::AbstractMeasure)
    FC, MA, MB = Core.Typeof(f_c), Core.Typeof(α), Core.Typeof(β)
    CombinedMeasure{FC,MA,MB}(f_c, α, β)
end

@inline function _generic_mcombine_impl_stage2(f_c, α::Dirac, β::Dirac)
    Dirac(f_c(α.x, β.x))
end


"""
    struct CombinedMeasure <: AbstractMeasure

Represents a combination of two measures.

User code should not create instances of `CombinedMeasure` directly, but
should call [`mcombine(f_c, α, β)`](@ref) instead.
"""
struct CombinedMeasure{FC,MA<:AbstractMeasure,MB<:AbstractMeasure} <: AbstractMeasure
    f_c::FC
    α::MA
    β::MB
end

@inline function preferred_stdmeasure(::Type{<:CombinedMeasure{<:Any,MA,MB}}) where {MA,MB}
    promote_stdmeasure(preferred_stdmeasure(MA), preferred_stdmeasure(MB))
end


@inline insupport(μ::CombinedMeasure, ab) = NoFastInsupport{typeof(μ)}()

@inline function mspace_flatsize(μ::CombinedMeasure{typeof(vcat)})
    _vcat_flatsize(mspace_flatsize(μ.α), mspace_flatsize(μ.β))
end

@inline _vcat_flatsize(a::SizeLike, b::SizeLike) = canonical_size((size2length(a) + size2length(b),))
@inline _vcat_flatsize(a::NoMSpaceElementSize, ::SizeLike) = a
@inline _vcat_flatsize(::SizeLike, b::NoMSpaceElementSize) = b
@inline _vcat_flatsize(a::NoMSpaceElementSize, ::NoMSpaceElementSize) = a

@inline getdof(μ::CombinedMeasure) = getdof(μ.α) + getdof(μ.β)
@inline fast_dof(μ::CombinedMeasure) = fast_dof(μ.α) + fast_dof(μ.β)

# Bypass `checked_arg`, would require splitting ab:
@inline checked_arg(::CombinedMeasure, ab) = ab

mdomain(μ::CombinedMeasure) = combinesets(μ.f_c, mdomain(μ.α), mdomain(μ.β))

rootmeasure(μ::CombinedMeasure) = mcombine(μ.f_c, rootmeasure(μ.α), rootmeasure(μ.β))

basemeasure(μ::CombinedMeasure) = mcombine(μ.f_c, basemeasure(μ.α), basemeasure(μ.β))

function logdensity_def(μ::CombinedMeasure, ab)
    # Use tpmeasure_split_combined to avoid duplicate calculation of transportmeasure(α):
    tpm_α, a, b = tpmeasure_split_combined(μ.f_c, μ.α, ab)
    return logdensity_def(tpm_α, a) + logdensity_def(μ.β, b)
end

# Density evaluation consumes the variate parts of both component measures
# in a single pass, using the with-rest protocol for value-dependent
# variate sizes:

logdensityof_impl(μ::CombinedMeasure, ab) = _combined_ld_impl(μ.f_c, μ, ab)

unsafe_logdensityof(μ::CombinedMeasure, ab) = logdensityof_impl(μ, ab)

function _combined_ld_impl(::typeof(tuple), μ::CombinedMeasure, ab::Tuple{Vararg{Any,2}})
    logdensityof(μ.α, ab[1]) + logdensityof(μ.β, ab[2])
end

function _combined_ld_impl(::Type{Pair}, μ::CombinedMeasure, ab::Pair)
    logdensityof(μ.α, ab.first) + logdensityof(μ.β, ab.second)
end

function _combined_ld_impl(::typeof(vcat), μ::CombinedMeasure, ab::AbstractVector)
    _point_result(_materialize(_combined_batched_ld(μ, ab, static(true))), μ)
end

function _combined_ld_impl(::typeof(merge), μ::CombinedMeasure, ab::NamedTuple)
    ℓ, _, x_rest = logdensityof_with_rest(μ, ab)
    isempty(x_rest) || _throw_stream_too_long()
    return ℓ
end

function _combined_ld_impl(f_c, μ::CombinedMeasure, ab)
    tpm_α, a, b = tpmeasure_split_combined(f_c, μ.α, ab)
    return logdensityof(tpm_α, a) + logdensityof(μ.β, b)
end

@inline mspace_ndims(::Type{<:CombinedMeasure{typeof(vcat)}}) = 1
@inline function fixed_stream_size(::Type{<:CombinedMeasure{<:Any,MA,MB}}) where {MA,MB}
    static(fixed_stream_size(MA) === static(true) && fixed_stream_size(MB) === static(true))
end

# Batches of vcat-combined variates are batches of streams: with fixed
# component sizes the whole batch is consumed in fused operations,
# otherwise stream by stream.
@inline function batched_logdensityof_impl(μ::CombinedMeasure{typeof(vcat)}, X::AbstractArray)
    _combined_batched_ld(μ, X, fixed_stream_size(μ))
end
function _combined_batched_ld(μ::CombinedMeasure, X::AbstractArray, ::True)
    ℓ, X_rest = batched_logdensityof_with_rest(μ, X, ())
    size(X_rest, 1) == 0 || _throw_stream_too_long()
    return ℓ
end
_combined_batched_ld(μ::CombinedMeasure, X::AbstractVector, ::False) = _combined_batched_ld(μ, X, static(true))
_combined_batched_ld(μ::CombinedMeasure, X::AbstractArray, ::False) = _streamwise_ld(logdensityof_impl, μ, X)

function batched_logdensityof_with_rest(μ::CombinedMeasure{typeof(vcat)}, X::AbstractArray, ::Tuple{})
    _combined_ld_with_rest(μ, X)
end
batched_logdensityof_with_rest(μ::CombinedMeasure{typeof(vcat)}, x::AbstractVector, ::Tuple{}) = _combined_ld_with_rest(μ, x)
function _combined_ld_with_rest(μ::CombinedMeasure, X::AbstractArray)
    ℓ_a, X2 = batched_logdensityof_with_rest(μ.α, X, ())
    ℓ_b, X_rest = batched_logdensityof_with_rest(μ.β, X2, ())
    return _lazy_add(ℓ_a, ℓ_b), X_rest
end

# Several variates per stream interleave the component parts, so the rows
# of each variate are split by the fixed component sizes:
function batched_logdensityof_with_rest(μ::CombinedMeasure{typeof(vcat)}, X::AbstractArray, sz::Dims)
    n_a, n_b = _fixed_stream_length(μ.α), _fixed_stream_length(μ.β)
    X_μ, X_rest = _batched_split(X, (n_a + n_b) * prod(sz))
    X_v = reshape(X_μ, (n_a + n_b, sz..., Base.tail(size(X_μ))...))
    X_a, X_b = _batched_split(X_v, n_a)
    ℓ_a, _ = batched_logdensityof_with_rest(μ.α, X_a, ())
    ℓ_b, _ = batched_logdensityof_with_rest(μ.β, X_b, ())
    return _lazy_add(ℓ_a, ℓ_b), X_rest
end

@inline _fixed_stream_length(μ) = _fixed_stream_length(μ, mspace_flatsize(μ))
@inline _fixed_stream_length(μ, sz::SizeLike) = dynamic(size2length(sz))
@noinline function _fixed_stream_length(μ, ::NoMSpaceElementSize)
    throw(ArgumentError("Consuming several variates per stream requires measures of type $(nameof(typeof(μ))) to have a known variate size"))
end

function logdensityof_with_rest(μ::CombinedMeasure{typeof(vcat)}, x::AbstractVector)
    ℓ_a, a, x2 = logdensityof_with_rest(μ.α, x)
    ℓ_b, b, x_rest = logdensityof_with_rest(μ.β, x2)
    x_μ, _ = _split_after(x, maybestatic_length(x) - maybestatic_length(x_rest))
    return ℓ_a + ℓ_b, x_μ, x_rest
end

function logdensityof_with_rest(μ::CombinedMeasure{typeof(merge)}, x::NamedTuple)
    ℓ_a, a, x2 = logdensityof_with_rest(μ.α, x)
    ℓ_b, b, x_rest = logdensityof_with_rest(μ.β, x2)
    return ℓ_a + ℓ_b, merge(a, b), x_rest
end


rand_impl(ctx::GenContext, μ::CombinedMeasure) = μ.f_c(rand_impl(ctx, μ.α), rand_impl(ctx, μ.β))

# Batches of vcat-combined measures are concatenated along the streams:
function batched_rand_impl(ctx::GenContext, μ::CombinedMeasure{typeof(vcat)}, sz::Dims)
    _combined_batched_rand(ctx, μ, sz, mspace_flatsize(μ.α), mspace_flatsize(μ.β))
end
function _combined_batched_rand(ctx::GenContext, μ::CombinedMeasure, sz::Dims, sz_a::SizeLike, sz_b::SizeLike)
    A = _as_stream_batch(batched_rand_impl(ctx, μ.α, sz), sz_a)
    B = _as_stream_batch(batched_rand_impl(ctx, μ.β, sz), sz_b)
    return vcat(A, B)
end
function _combined_batched_rand(ctx::GenContext, μ::CombinedMeasure, sz::Dims, ::Any, ::Any)
    _batched_rand_pointwise(ctx, μ, sz)
end


# Transport consumes the variate parts of both component measures in a
# single pass, analogous to density evaluation:

transport_to_std(::Type{S}, μ::CombinedMeasure, ab) where {S<:StdMeasure} = _combined_to_std(S, μ.f_c, μ, ab)

function _combined_to_std(::Type{S}, f_c, μ::CombinedMeasure, ab) where {S}
    tpm_α, a, b = tpmeasure_split_combined(f_c, μ.α, ab)
    vcat(_as_stdstream(transport_to_std(S, tpm_α, a)), _as_stdstream(transport_to_std(S, μ.β, b)))
end

function _combined_to_std(::Type{S}, ::Union{typeof(vcat),typeof(merge)}, μ::CombinedMeasure, ab) where {S}
    z, _, x_rest = transport_to_std_with_rest(S, μ, ab)
    if !isempty(x_rest)
        throw(ArgumentError("Variate too long during transport of a combined measure"))
    end
    return z
end

function transport_to_std_with_rest(::Type{S}, μ::CombinedMeasure{typeof(vcat)}, x::AbstractVector) where {S<:StdMeasure}
    z_a, _, x2 = transport_to_std_with_rest(S, μ.α, x)
    z_b, _, x_rest = transport_to_std_with_rest(S, μ.β, x2)
    x_μ, _ = _split_after(x, maybestatic_length(x) - maybestatic_length(x_rest))
    return vcat(z_a, z_b), x_μ, x_rest
end

function transport_to_std_with_rest(::Type{S}, μ::CombinedMeasure{typeof(merge)}, x::NamedTuple) where {S<:StdMeasure}
    z_a, a, x2 = transport_to_std_with_rest(S, μ.α, x)
    z_b, b, x_rest = transport_to_std_with_rest(S, μ.β, x2)
    return vcat(z_a, z_b), merge(a, b), x_rest
end

function transport_from_std_with_rest(::Type{S}, μ::CombinedMeasure, z::AbstractVector) where {S<:StdMeasure}
    a, z2 = transport_from_std_with_rest(S, μ.α, z)
    b, z_rest = transport_from_std_with_rest(S, μ.β, z2)
    return μ.f_c(a, b), z_rest
end


# Batched transport consumes the variate parts of both component measures
# along batches of streams:

function batched_transport_to_std(::Type{S}, μ::CombinedMeasure{typeof(vcat)}, X::AbstractArray) where {S<:StdMeasure}
    Z, _, X_rest = batched_transport_to_std_with_rest(S, μ, X)
    if size(X_rest, 1) != 0
        throw(ArgumentError("Variate streams too long during batched transport of a combined measure"))
    end
    return Z
end

function batched_transport_to_std_with_rest(::Type{S}, μ::CombinedMeasure{typeof(vcat)}, X::AbstractArray) where {S<:StdMeasure}
    Z_a, _, X2 = batched_transport_to_std_with_rest(S, μ.α, X)
    Z_b, _, X_rest = batched_transport_to_std_with_rest(S, μ.β, X2)
    X_μ, _ = _batched_split(X, size(X, 1) - size(X_rest, 1))
    return vcat(Z_a, Z_b), X_μ, X_rest
end

function batched_transport_from_std(::Type{S}, μ::CombinedMeasure{typeof(vcat)}, Z::AbstractArray) where {S<:StdMeasure}
    X, Z_rest = batched_transport_from_std_with_rest(S, μ, Z)
    if size(Z_rest, 1) != 0
        throw(ArgumentError("Length of standard variates doesn't match degrees of freedom of a combined measure"))
    end
    return X
end

function batched_transport_from_std_with_rest(::Type{S}, μ::CombinedMeasure{typeof(vcat)}, Z::AbstractArray) where {S<:StdMeasure}
    A, Z2 = batched_transport_from_std_with_rest(S, μ.α, Z)
    B, Z_rest = batched_transport_from_std_with_rest(S, μ.β, Z2)
    X = vcat(_as_stream_batch(A, mspace_flatsize(μ.α)), _as_stream_batch(B, mspace_flatsize(μ.β)))
    return X, Z_rest
end
