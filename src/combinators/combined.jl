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

_split_variate_byvalue(::typeof(vcat), ::Real, ab::AbstractVector) =
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

_mcombine_product_shortcut(::typeof(vcat), ma::AbstractVector, mb::AbstractVector, α, β) =
    productmeasure(vcat(ma, mb))
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


@inline insupport(μ::CombinedMeasure, ab) = NoFastInsupport{typeof(μ)}()

@inline getdof(μ::CombinedMeasure) = getdof(μ.α) + getdof(μ.β)
@inline fast_dof(μ::CombinedMeasure) = fast_dof(μ.α) + fast_dof(μ.β)

# Bypass `checked_arg`, would require splitting ab:
@inline checked_arg(::CombinedMeasure, ab) = ab

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


function transport_to_mvstd(ν_inner::StdMeasure, μ::CombinedMeasure, ab)
    tpm_α, a, b = tpmeasure_split_combined(μ.f_c, μ.α, ab)
    y1 = transport_to_mvstd(ν_inner, tpm_α, a)
    y2 = transport_to_mvstd(ν_inner, μ.β, b)
    return vcat(y1, y2)
end


function transport_from_mvstd_with_rest(ν::CombinedMeasure, μ_inner::StdMeasure, x)
    a, x2 = transport_from_mvstd_with_rest(ν.α, μ_inner, x)
    b, x_rest = transport_from_mvstd_with_rest(ν.β, μ_inner, x2)
    return ν.f_c(a, b), x_rest
end
