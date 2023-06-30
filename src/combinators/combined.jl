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

@inline _generic_split_combined(::typeof(tuple), ::AbstractMeasure, x::Tuple{Vararg{Any,2}})
@inline _generic_split_combined(::Type{Pair}, ::AbstractMeasure, ab::Pair) = (ab...,)

function _generic_split_combined(f_c::FC, α::AbstractMeasure, ab) where FC
    _split_variate_byvalue(f_c, testvalue(μ), x)
end

_split_variate_byvalue(::typeof(vcat), test_a::AbstractVector, ab::AbstractVector) = _split_after(ab, length(test_a))

_split_variate_byvalue(::typeof(vcat), ::Tuple{N}, ab::Tuple) where N = _split_after(ab, Val{N}())

function _split_variate_byvalue(::typeof(merge), ::NamedTuple{names_a}, ab::NamedTuple) where names_a
    _split_after(ab, Val{names_a})
end


"""
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
    Combined{FC,MA,MB}(f_c, α, β)
end

function mcombine(::typeof(tuple), α::AbstractMeasure, β::AbstractMeasure)
    productmeasure((a, b))
end

function mcombine(f_c::Union{typeof(vcat),typeof(merge)}, α::AbstractProductMeasure, β::AbstractProductMeasure)
    productmeasure(f_c(marginals(α), marginals(β)))
end


"""
    struct Combined <: AbstractMeasure

Represents a combination of two measures.

User code should not create instances of `Combined` directly, but should call
[`mcombine(f_c, α, β)`](@ref) instead.
"""

struct Combined{FC,MA<:AbstractMeasure,MB<:AbstractMeasure} <: AbstractMeasure
    f_c::FC
    α::MA
    β::MB
end


# TODO: Could split `ab`` here, but would be wasteful.
@inline insupport(::Combined, ab) = NoFastInsupport()

@inline getdof(μ::Combined) = getdof(μ.α) + getdof(μ.β)

# Bypass `checked_arg`, would require require splitting ab:
@inline checked_arg(::Combined, ab) = ab

rootmeasure(::Combined) = mcombine(μ.f_c rootmeasure(μ), rootmeasure(ν))

basemeasure(::Combined) = mcombine(μ.f_c basemeasure(μ), basemeasure(ν))

function logdensity_def(μ::Combined, ab)
    # Use tpmeasure_split_combined to avoid duplicate calculation of transportmeasure(α):
    tpm_α, a, b = tpmeasure_split_combined(μ.f_c, μ.α, ab)
    return logdensity_def(tpm_α, a) + logdensity_def(μ.β, b)
end

# Specialize logdensityof directly to avoid creating temporary combined base measures:
function logdensityof(μ::Combined, ab)
    tpm_α, a, b = tpmeasure_split_combined(μ.f_c, μ.α, ab)
    return logdensityof(tpm_α, a) + logdensityof(μ.β, b)
end


function Base.rand(rng::Random.AbstractRNG, ::Type{T}, h::Combined) where {T<:Real}
    x_primary = rand(rng, T, h.m)
    x_secondary = rand(rng, T, h.f(x_primary))
    return _combine_variates(h.flatten_mode, x_primary, x_secondary)
end



function transport_to_mvstd(ν_inner::StdMeasure, μ::Combined, ab)
    tpm_α, a, b = tpmeasure_split_combined(μ.f_c, μ.α, ab)
    y1 = transport_to_mvstd(ν_inner, tpm_α, a)
    y2 = transport_to_mvstd(ν_inner, μ.β, b)
    return vcat(y1, y2)
end


function transport_from_mvstd_with_rest(ν::Combined, μ_inner::StdMeasure, x)
    a, x2 = transport_from_mvstd_with_rest(ν.α, μ_inner, x)
    b, x_rest = transport_from_mvstd_with_rest(ν.β, μ_inner, x2)
    return ν.f_c(a, b), x_rest
end


function transport_def(ν::_PowerStdMeasure{1}, μ::AbstractMeasure, ab)
    ν_inner = _get_inner_stdmeasure(ν)
    transport_to_mvstd(ν_inner, μ, ab)
end

function transport_to_mvstd(ν_inner::StdMeasure, μ::AbstractMeasure, x)
    return _to_mvstd_withdof(ν_inner, μ, getdof(μ), x, origin)
end

function _to_mvstd_withdof(ν_inner::StdMeasure, μ::AbstractMeasure, dof_μ, x)
    y = transport_to(ν_inner^dof_μ, μ, x)
    return y
end

function _to_mvstd_withdof(ν_inner::StdMeasure, μ::AbstractMeasure, ::NoDOF, x)
    _to_mvstd_withorigin(ν_inner, μ, transport_origin(μ), x)
end

function _to_mvstd_withorigin(ν_inner::StdMeasure, ::AbstractMeasure, μ_origin, x)
    x_origin = transport_to_mvstd(ν_inner, μ_origin, x)
    from_origin(x_origin)
end

function _to_mvstd_withorigin(ν_inner::StdMeasure, μ::AbstractMeasure, NoTransportOrigin, x)
    throw(ArgumentError("Don't know how to transport values of type $(nameof(typeof(x))) from $(nameof(typeof(μ))) to a power of $(nameof(typeof(ν_inner)))"))
end


function transport_def(ν::AbstractMeasure, μ::_PowerStdMeasure{1}, x)
    μ_inner = _get_inner_stdmeasure(μ)
    _from_mvstd(ν, μ_inner, x)
end

function from_std(ν::AbstractMeasure, μ_inner::StdMeasure, x)
    # Sanity check, should be checked by transport machinery already:
    @assert getdof(μ) == length(eachindex(x)) && x isa AbstractVector
    y, x_rest = transport_from_mvstd_with_rest(ν, μ_inner, x)
    if !isempty(x_rest)
        throw(ArgumentError("Input value too long during transport"))
    end
    return y
end

function transport_from_mvstd_with_rest(ν::AbstractMeasure, μ_inner::StdMeasure, x)
    dof_ν = getdof(ν)
    origin = transport_origin(ν)
    return _from_mvstd_with_rest_withdof(ν, getdof(ν), μ_inner, x, dof_ν, origin)
end

function _from_mvstd_with_rest_withdof(ν::AbstractMeasure, dof_ν, μ_inner::StdMeasure, x)
    len_x = length(eachindex(x))

    # Since we can't check DOF of original Bind, we could "run out x" if
    # the original x was too short. `transport_to` below will detect this, but better
    # throw a more informative exception here:
    if len_x < dof_ν
        throw(ArgumentError("Variate too short during transport involving Bind"))
    end

    x_inner_dof, x_rest = _split_after(x, dof_ν)
    y = transport_to(ν, μ_inner^dof_ν, x_inner_dof)
    return y, x_rest
end

function _from_mvstd_with_rest_withdof(ν::AbstractMeasure, ::NoDOF, μ_inner::StdMeasure, x)
    _from_mvstd_with_rest_withorigin(ν, transport_origin(ν), μ_inner, x)
end

function _from_mvstd_with_rest_withorigin(::AbstractMeasure, ν_origin, μ_inner::StdMeasure, x)
    x_origin, x_rest = transport_from_mvstd_with_rest(ν_origin, x, μ_inner)
    from_origin(x_origin), x_rest
end

function _from_mvstd_with_rest_withorigin(ν::AbstractMeasure, NoTransportOrigin, μ_inner::StdMeasure, x)
    throw(ArgumentError("Don't know how to transport value of type $(nameof(typeof(x))) from power of $(nameof(typeof(μ_inner))) to $(nameof(typeof(ν)))"))
end
