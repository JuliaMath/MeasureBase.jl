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

function _generic_split_combined(f_c::FC, α::AbstractMeasure, ab) where FC
    _split_variate_byvalue(f_c, testvalue(μ), ab)
end

_split_variate_byvalue(::typeof(vcat), test_a::AbstractVector, ab::AbstractVector) = _split_after(ab, length(test_a))

_split_variate_byvalue(::typeof(vcat), ::Tuple{N}, ab::Tuple) where N = _split_after(ab, Val{N}())

function _split_variate_byvalue(::typeof(merge), ::NamedTuple{names_a}, ab::NamedTuple) where names_a
    _split_after(ab, Val{names_a})
end


@doc raw"""
    mcombine(f_c, α::AbstractMeasure, β::AbstractMeasure)

Combines two measures `α` and `β` to a combined measure via a point combination
function `f_c`.

`f_c` must combine a given point `a` from the space of measure `α` with a
given point `b` from the space of measure `β` to a single value
`ab = f_c(a, b)` in the space of the combined measure
`μ = mcombine(f_c, α, β)`.

The combined measure has the mathethematical interpretation (on
sets $$A$$ and $$B$$)

```math
\mu(f_c(A, B)) = \alpha(A)\, \beta(B)
```
"""
function mcombine end
export mcombine

function mcombine(f_c, α::AbstractMeasure, β::AbstractMeasure)
    FC, MA, MB = Core.Typeof(f_c), Core.Typeof(α), Core.Typeof(β)
    CombinedMeasure{FC,MA,MB}(f_c, α, β)
end

function mcombine(::typeof(tuple), α::AbstractMeasure, β::AbstractMeasure)
    productmeasure((a, b))
end

function mcombine(f_c::Union{typeof(vcat),typeof(merge)}, α::AbstractProductMeasure, β::AbstractProductMeasure)
    productmeasure(f_c(marginals(α), marginals(β)))
end


"""
    struct CombinedMeasure <: AbstractMeasure

Represents a combination of two measures.

User code should not create instances of `CombinedMeasure` directly, but should call
[`mcombine(f_c, α, β)`](@ref) instead.
"""

struct CombinedMeasure{FC,MA<:AbstractMeasure,MB<:AbstractMeasure} <: AbstractMeasure
    f_c::FC
    α::MA
    β::MB
end


@inline insupport(μ::CombinedMeasure, ab) = NoFastInsupport{typeof(μ)}()

@inline getdof(μ::CombinedMeasure) = getdof(μ.α) + getdof(μ.β)
@inline fast_dof(μ::CombinedMeasure) = fast_dof(μ.α) + fast_dof(μ.β)

# Bypass `checked_arg`, would require require splitting ab:
@inline checked_arg(::CombinedMeasure, ab) = ab

rootmeasure(::CombinedMeasure) = mcombine(μ.f_c, rootmeasure(μ), rootmeasure(ν))

basemeasure(::CombinedMeasure) = mcombine(μ.f_c, basemeasure(μ), basemeasure(ν))

function logdensity_def(μ::CombinedMeasure, ab)
    # Use tpmeasure_split_combined to avoid duplicate calculation of transportmeasure(α):
    tpm_α, a, b = tpmeasure_split_combined(μ.f_c, μ.α, ab)
    return logdensity_def(tpm_α, a) + logdensity_def(μ.β, b)
end

# Specialize logdensityof directly to avoid creating temporary combined base measures:
function logdensityof(μ::CombinedMeasure, ab)
    tpm_α, a, b = tpmeasure_split_combined(μ.f_c, μ.α, ab)
    return logdensityof(tpm_α, a) + logdensityof(μ.β, b)
end


function Base.rand(rng::Random.AbstractRNG, ::Type{T}, h::CombinedMeasure) where {T<:Real}
    x_primary = rand(rng, T, h.m)
    x_secondary = rand(rng, T, h.f(x_primary))
    return _combine_variates(h.flatten_mode, x_primary, x_secondary)
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
