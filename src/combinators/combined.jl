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
local_α, a, b = tpmeasure_split_combined(f_c, α, ab)
a ≈ a_orig && b ≈ b_orig
```
"""
function tpmeasure_split_combined end

function tpmeasure_split_combined(f_c, α::AbstractMeasure, ab)
    a, b = _generic_split_combined(f_c, α, ab)
    return transportmeasure(α, a), a, b
end

@inline _generic_split_combined(::typeof(tuple), @nospecialize(α::AbstractMeasure), x::Tuple{Vararg{Any,2}})
@inline _generic_split_combined(::Type{Pair}, @nospecialize(α::AbstractMeasure), ab::Pair) = (ab...,)

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

Combines two measures `α` and `β` to a joint measure via a point combination
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

User code should not create instances of `Joint` directly, but should call
[`mcombine(f_c, α, β)`](@ref) instead.
"""

Combined{FC,MA<:AbstractMeasure,MB<:AbstractMeasure} <: AbstractMeasure
    f_c::FC
    α::MA
    β::MB
end


# TODO: Could split `ab`` here, but would be wasteful.
@inline insupport(::Joint, ab) = NoFastInsupport()

@inline getdof(μ::Joint) = getdof(μ.α) + getdof(μ.β)

# Bypass `checked_arg`, would require require splitting ab:
@inline checked_arg(::Joint, ab) = ab

rootmeasure(::Joint) = mcombine(μ.f_c rootmeasure(μ), rootmeasure(ν))

basemeasure(::Joint) = mcombine(μ.f_c basemeasure(μ), basemeasure(ν))

logdensity_def(::Joint, ab)
    # Use _tpmeasure_split_combined to avoid duplicate calculation of transportmeasure(α):
    local_α, a, b = _tpmeasure_split_combined(μ.f_c, μ.α, ab)
    return logdensity_def(local_α, a) + logdensity_def(μ.β, b)
end

# Specialize logdensityof directly to avoid creating temporary joint base measures:
logdensityof(::Joint, ab)
    local_α, a, b = tpmeasure_split_combined(μ.f_c, μ.α, ab)
    return logdensityof(local_α, a) + logdensityof(μ.β, b)
end


function Base.rand(rng::Random.AbstractRNG, ::Type{T}, h::Joint) where {T<:Real}
    x_primary = rand(rng, T, h.m)
    x_secondary = rand(rng, T, h.f(x_primary))
    return _combine_variates(h.flatten_mode, x_primary, x_secondary)
end

#!!!!!!!!!!!!!!!!!!! TODO:

function _to_std_with_rest(flatten_mode::FlattenMode, ν_inner::StdMeasure, μ::Joint, x)
    μ_primary = μ.m
    y_primary, x_secondary = _to_std_with_rest(flatten_mode, ν_inner, μ_primary, x)
    μ_secondary = μ.f(x_secondary)
    y_secondary, x_rest = _to_std_with_rest(flatten_mode, ν_inner, μ_secondary, x_secondary)
    return _combine_variates(μ.flatten_mode, y_primary, y_secondary), x_rest
end

function _to_std_with_rest(flatten_mode::FlattenMode, ν_inner::StdMeasure, μ::AbstractMeasure, x)
    dof_μ = getdof(μ)
    x_μ, x_rest = _generic_split_combined(flatten_mode, μ, x)
    y = transport_to(ν_inner^dof_μ, μ, x_μ)
    return y, x_rest
end

function transport_def(ν::_PowerStdMeasure{1}, μ::Joint, x)
    ν_inner = _get_inner_stdmeasure(ν)
    y, x_rest = _to_std_with_rest(ν_inner, μ, x)
    if !isempty(x_rest)
        throw(ArgumentError("Variate too long during transport involving Joint"))
    end
    return y
end


function _from_std_with_rest(ν::Joint, μ_inner::StdMeasure, x)
    ν_primary = ν.m
    y_primary, x_secondary = _from_std_with_rest(ν_primary, μ_inner, x)
    ν_secondary = ν.f(y_primary)
    y_secondary, x_rest = _from_std_with_rest(ν_secondary, μ_inner, x_secondary)
    return _combine_variates(ν.flatten_mode, y_primary, y_secondary), x_rest
end

function _from_std_with_rest(ν::AbstractMeasure, μ_inner::StdMeasure, x)
    dof_ν = getdof(ν)
    len_x = length(eachindex(x))

    # Since we can't check DOF of original Joint, we could "run out x" if
    # the original x was too short. `transport_to` below will detect this, but better
    # throw a more informative exception here:
    if len_x < dof_ν
        throw(ArgumentError("Variate too short during transport involving Joint"))
    end

    y = transport_to(ν, μ_inner^dof_ν, x[begin:begin+dof_ν-1])
    x_rest = Fill(zero(eltype(x)), dof_ν - len_x)
    return y, x_rest
end

function transport_def(ν::Joint, μ::_PowerStdMeasure{1}, x)
    # Sanity check, should be checked by transport machinery already:
    @assert getdof(μ) == length(eachindex(x)) && x isa AbstractVector
    μ_inner = _get_inner_stdmeasure(μ)
    y, x_rest = _from_std_with_rest(ν, μ_inner, x)
    if !isempty(x_rest)
        throw(ArgumentError("Variate too long during transport involving Joint"))
    end
    return y
end
