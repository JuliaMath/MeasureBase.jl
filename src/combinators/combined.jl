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

function _combined_ld_impl(::Union{typeof(vcat),typeof(merge)}, μ::CombinedMeasure, ab)
    ℓ, x_μ, x_rest = logdensityof_with_rest(μ, ab)
    if !isempty(x_rest)
        throw(
            ArgumentError(
                "Variate too long during density evaluation of a combined measure",
            ),
        )
    end
    return ℓ
end

function _combined_ld_impl(f_c, μ::CombinedMeasure, ab)
    tpm_α, a, b = tpmeasure_split_combined(f_c, μ.α, ab)
    return logdensityof(tpm_α, a) + logdensityof(μ.β, b)
end

function batched_logdensityof_impl(μ::CombinedMeasure{typeof(vcat)}, A::AbstractArray)
    ℓ, _, A_rest = batched_logdensityof_with_rest(μ, A)
    if size(A_rest, 1) != 0
        throw(ArgumentError("Variate streams too long during batched density evaluation of a combined measure"))
    end
    return ℓ
end

function batched_logdensityof_with_rest(μ::CombinedMeasure{typeof(vcat)}, A::AbstractArray)
    ℓ_a, _, A2 = batched_logdensityof_with_rest(μ.α, A)
    ℓ_b, _, A_rest = batched_logdensityof_with_rest(μ.β, A2)
    n_μ = size(A, 1) - size(A_rest, 1)
    A_μ = view(A, 1:n_μ, Base.tail(axes(A))...)
    return _lazy_add(ℓ_a, ℓ_b), A_μ, A_rest
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


function Base.rand(rng::Random.AbstractRNG, ::Type{T}, μ::CombinedMeasure) where {T<:Real}
    a = rand(rng, T, μ.α)
    b = rand(rng, T, μ.β)
    return μ.f_c(a, b)
end

function Base.rand(rng::Random.AbstractRNG, μ::CombinedMeasure)
    a = rand(rng, μ.α)
    b = rand(rng, μ.β)
    return μ.f_c(a, b)
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
