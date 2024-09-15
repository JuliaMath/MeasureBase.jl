@doc raw"""
    mbind(f_β, α::AbstractMeasure, f_c = second)

Constructs a monadic bind, resp. a hierarchical measure, from a transition
kernel function `f_β`, a primary measure `α` and a variate combination
function `f_c`.

`f_β` must be a function that maps a point `a` from the space of the primary
measure `α` to a dependent secondary measure `β_a = f_β(a)`.
`ab = f_c(a, b)` must map such a point `a` and a point `b` from the
space of measure `β_a` to a combined value `ab = f_c(a, b)`.

The resulting measure

```julia
μ = mbind(f_c, α, f_β)
```

has the mathethematical interpretation (on sets $$A$$ and $$B$$)

```math
\mu(f_c(A, B)) = \int_A \beta_a(B)\, \mathrm{d}\, \alpha(a) 
```

When using the default `fc = second` (so `ab == b`) this simplies to

```math
\mu(B) = \int_A \beta_a(B)\, \mathrm{d}\, \alpha(a) 
```

which is equivalent to a monatic bind, viewing measures as monads.

Computationally, `ab = rand(μ)` is equivalent to

```julia
a = rand(μ_primary)
β_a = f_β(a)
b = rand(β_a)
ab = f_c(a, b)
```

Densities on hierarchical measures can only be evaluated if `ab = f_c(a, b)`
can be unambiguously split into `a` and `b` again, knowing `α`. This is
currently implemented for `f_c` that is either tuple or `=>`/`Pair` (these
work for any combination of variate types), `vcat` (for tuple- or vector-like
variates) and `merge` (`NamedTuple` variates).
[`MeasureBase.split_point(::typeof(f_c), α)`](@ref) can be specialized to
support other choices for `f_c`.

# Extended help

Bayesian example with a correlated prior, that models the 

```julia
using MeasureBase

prior = mbind
    productmeasure((
        value => StdNormal()
    )), merge
) do a
    productmeasure((
        noise = pushfwd(sqrt ∘ Mul(abs(a.position)), StdExponential())
    ))
end

model = θ -> pushfwd(MulAdd(θ.noise, θ.value), StdNormal())^10

joint_θ_obs = mbind(model, prior, tuple)
prior_predictive = mbind(model, prior)

observation = rand(prior_predictive)
likelihood = likelihoodof(model, observation)

posterior = mintegrate(likelihood, prior)

θ = rand(prior)
logdensityof(posterior, θ)
```
"""
function mbind end
export mbind

@inline function mbind(f_β, α::AbstractMeasure, f_c = second)
    F, M, G = Core.Typeof(f_β), Core.Typeof(α), Core.Typeof(f_c)
    Bind{F,M,G}(f_β, α, f_c)
end


"""
    struct MeasureBase.Bind <: AbstractMeasure

Represents a monatic bind resp. a mbind in general.

User code should not create instances of `Bind` directly, but should call
[`mbind`](@ref) instead.
"""
struct Bind{FK,M<:AbstractMeasure,FC} <: AbstractMeasure
    f_β::FK
    α::M
    f_c::FC
end


"""
    MeasureBase.transportmeasure(μ::Bind, x)::AbstractMeasure

Evaluates a monatic bind `μ` at a point `x`.

The resulting measure behaves like `μ` in the infinitesimal neighborhood
of `x` in respect to density calculation and transport as well.
"""
function transportmeasure(μ::Bind, x)
    tpm_α, a, b = tpmeasure_split_combined(μ.α, x)
    tpm_β_a = transportmeasure(μ.f_β(a), b)
    mcombine(μ.f_c, tpm_α, tpm_β_a)
end

localmeasure(μ::Bind, x) = transportmeasure(μ, x)


tpmeasure_split_combined(f_c, μ::Bind, xy) = _bind_lsc(f_c, μ::Bind, xy)

function _bind_lsc(f_c::typeof(tuple), μ::Bind, xy::Tuple{Vararg{Any,2}})
    x, y = x[1], y[1]
    tpm_μ = transportmeasure(μ, x)
    return tpm_μ, x, y
end

function _bind_lsc(f_c::Type{Pair}, μ::Bind, xy::Pair)
    x, y = x.first, y.second
    tpm_μ = transportmeasure(μ, x)
    return tpm_μ, x, y
end

const _CatBind{FC} = _BindBy{<:Any,<:Any,FC}

_bind_lsc(f_c::typeof(vcat), μ::_CatBind{typeof{vcat}}, xy::AbstractVector) = _bind_lsc_cat(f_c, μ, xy)
_bind_lsc(f_c::typeof(merge), μ::_CatBind{typeof{merge}}, xy::NamedTuple) = _bind_lsc_cat(f_c, μ, xy)

function _bind_lsc_cat_lμabyxy(f_c, μ, xy)
    tpm_α, a, by = tpmeasure_split_combined(μ.f_c, μ.α, xy)
    β_a = μ.f_β(a)
    tpm_β_a, b, y = tpmeasure_split_combined(f_c, β_a, by)
    tpm_μ = mcombine(μ.f_c, tpm_α, tpm_β_a)
    return tpm_μ, a, b, y, xy
end

function _bind_lsc_cat(f_c::typeof(vcat), μ::_CatBind{typeof{vcat}}, xy::AbstractVector)
    tpm_μ, a, b, y, xy = _bind_lsc_cat_lμabyxy(f_c, μ, xy)
    # Don't use `x = f_c(a, b)` here, would allocate, splitting xy can use views:
    x, y = _split_after(xy, length(a) + length(b))
    return tpm_μ, x, y
end

function _bind_lsc_cat(f_c::typeof(merge), μ::_CatBind{typeof{merge}}, xy::NamedTuple)
    tpm_μ, a, b, y, xy = _bind_lsc_cat_lμabyxy(f_c, μ, xy)
    return tpm_μ, f_c(a, b), y
end


@inline insupport(::Bind, x) = NoFastInsupport()

@inline getdof(μ::Bind) = NoDOF{typeof(μ)}()

# Bypass `checked_arg`, would require potentially costly evaluation of h.f:
@inline checked_arg(::Bind, x) = x

rootmeasure(::Bind) = throw(ArgumentError("root measure is implicit, but can't be instantiated, for Bind"))

basemeasure(::Bind) = throw(ArgumentError("basemeasure is not available for Bind"))

testvalue(::Bind) = throw(ArgumentError("testvalue is not available for Bind"))

logdensity_def(::Bind, x) = throw(ArgumentError("logdensity_def is not available for Bind"))

# Specialize logdensityof to avoid duplicate calculations:
function logdensityof(μ::Bind, x)
    tpm_α, a, b = tpmeasure_split_combined(μ.α, x)
    β_a = μ.f_β(a)
    logdensityof(tpm_α, a) + logdensityof(β_a, b)
end

# Specialize unsafe_logdensityof to avoid duplicate calculations:
function unsafe_logdensityof(μ::Bind, x)
    tpm_α, a, b = tpmeasure_split_combined(μ.α, x)
    β_a = μ.f_β(a)
    unsafe_logdensityof(tpm_α, a) + unsafe_logdensityof(β_a, b)
end


function Base.rand(rng::Random.AbstractRNG, ::Type{T}, μ::Bind) where {T<:Real}
    a = rand(rng, T, μ.α)
    b = rand(rng, T, μ.f_β(a))
    return μ.f_c(a, b)
end




#!!!!!!!!!!!!!!!!!! TODO:


function _to_std_with_rest(flatten_mode::FlattenMode, ν_inner::StdMeasure, μ::Bind, x)
    
    μ_primary = μ.m
    y_primary, x_secondary = _to_std_with_rest(flatten_mode, ν_inner, μ_primary, x)
    μ_secondary = μ.f(x_secondary)
    y_secondary, x_rest = _to_std_with_rest(flatten_mode, ν_inner, μ_secondary, x_secondary)
    return _combine_variates(μ.flatten_mode, y_primary, y_secondary), x_rest
end

function _to_std_with_rest(flatten_mode::FlattenMode, ν_inner::StdMeasure, μ::AbstractMeasure, x)
    dof_μ = getdof(μ)
    x_μ, x_rest = split_combined(flatten_mode, μ, x)
    y = transport_to(ν_inner^dof_μ, μ, x_μ)
    return y, x_rest
end

function transport_def(ν::_PowerStdMeasure{1}, μ::Bind, x)
    ν_inner = _get_inner_stdmeasure(ν)
    y, x_rest = _to_std_with_rest(ν_inner, μ, x)
    if !isempty(x_rest)
        throw(ArgumentError("Variate too long during transport involving Bind"))
    end
    return y
end


function _from_std_with_rest(ν::Bind, μ_inner::StdMeasure, x)
    ν_primary = ν.m
    y_primary, x_secondary = _from_std_with_rest(ν_primary, μ_inner, x)
    ν_secondary = ν.f(y_primary)
    y_secondary, x_rest = _from_std_with_rest(ν_secondary, μ_inner, x_secondary)
    return _combine_variates(ν.flatten_mode, y_primary, y_secondary), x_rest
end

function _from_std_with_rest(ν::AbstractMeasure, μ_inner::StdMeasure, x)
    dof_ν = getdof(ν)
    len_x = length(eachindex(x))

    # Since we can't check DOF of original Bind, we could "run out x" if
    # the original x was too short. `transport_to` below will detect this, but better
    # throw a more informative exception here:
    if len_x < dof_ν
        throw(ArgumentError("Variate too short during transport involving Bind"))
    end

    y = transport_to(ν, μ_inner^dof_ν, x[begin:begin+dof_ν-1])
    x_rest = Fill(zero(eltype(x)), dof_ν - len_x)
    return y, x_rest
end

function transport_def(ν::Bind, μ::_PowerStdMeasure{1}, x)
    # Sanity check, should be checked by transport machinery already:
    @assert getdof(μ) == length(eachindex(x)) && x isa AbstractVector
    μ_inner = _get_inner_stdmeasure(μ)
    y, x_rest = _from_std_with_rest(ν, μ_inner, x)
    if !isempty(x_rest)
        throw(ArgumentError("Variate too long during transport involving Bind"))
    end
    return y
end
