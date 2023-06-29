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


_local_productmeasure(fc, α, \) = productmeasure(fc(marginals(μ1), marginals(μ2)))

# TODO: _local_productmeasure(::AutoFlatten, μ1, μ2) = productmeasure(μ1, μ2)
# Needs a FlatProductMeasure type.

function _local_split_combined(f_c, μ::Bind, x)
    local_α, a, b_withrest = _local_split_combined(μ.α, x)
    β_a = μ.f_β(a)
    local_β_a, b, rest = _localmeasure_with_rest(β_a, b_withrest)
    local_μ = productmeasure(μ.f_c(marginals(local_α), marginals(local_β_a)))
    local_ab, rest = split_combined(μ.f_c, local_μ, ab)
    return local_μ, local_ab, rest
end

function _local_split_combined(f_c, μ::AbstractMeasure, x_withrest)
    x, rest = split_combined(f_c, μ, x_withrest)
    return localmeasure(μ, x), x, rest
end


@inline insupport(::Bind, x) = NoFastInsupport()

@inline getdof(μ::Bind) = NoDOF{typeof(μ)}()

# Bypass `checked_arg`, would require potentially costly evaluation of h.f:
@inline checked_arg(::Bind, x) = x

rootmeasure(::Bind) = throw(ArgumentError("root measure is implicit, but can't be instantiated, for Bind"))

basemeasure(::Bind) = throw(ArgumentError("basemeasure is not available for Bind"))

logdensity_def(::Bind, x) = throw(ArgumentError("logdensity_def is not available for Bind"))


# # TODO: Default implementation of unsafe_logdensityof is a bit inefficient
# # for AutoFlatten, since variate will be split in `localmeasure` and then
# # split again in log-density evaluation. Maybe add something like
# function unsafe_logdensityof(h::Bind, x)
#     local_primary, local_secondary, x_primary, x_secondary = ...
#     # Need to call full logdensityof for h_secondary since x_secondary hasn't
#     # been checked yet:
#     unsafe_logdensityof(local_primary, x_primary) + logdensityof(local_secondary, x_secondary)
# end


function Base.rand(rng::Random.AbstractRNG, ::Type{T}, h::Bind) where {T<:Real}
    x_primary = rand(rng, T, h.m)
    x_secondary = rand(rng, T, h.f(x_primary))
    return _combine_variates(h.flatten_mode, x_primary, x_secondary)
end


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
