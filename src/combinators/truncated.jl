export TruncatedMeasure

"""
    struct TruncatedMeasure <: AbstractMeasure

A measure with scalar variates restricted to an interval, see
[`truncate`](@ref Base.truncate(::AbstractMeasure, ::AbstractInterval)).
"""
struct TruncatedMeasure{M,T,W} <: AbstractMeasure
    parent::M
    lo::T
    hi::T
    mass::W
end

"""
    truncate(μ::AbstractMeasure, a..b)

The measure `μ` restricted to the interval `a..b`, a
[`TruncatedMeasure`](@ref).

`μ` must have scalar variates, no atoms and a Stieltjes measure function
(see [`smf`](@ref)). The mass of the truncation is the mass of `μ` on the
interval, `smf(μ, a, b)`, and `normalize(truncate(μ, a..b))` is the
law of `μ` conditioned on the interval. Open and closed ends are treated
alike (they differ by null sets), the interval may be half-infinite. For
measures `μ` symmetric about zero, `normalize(truncate(μ, 0..Inf))` has the
law of [`Half(μ)`](@ref Half).

Densities are those of `μ` inside the interval and `-Inf` outside. Variates
are generated and transported via the interval masses `smf(μ, a, x)` of `μ`
and their inverse `invsmf(μ, a, p)` (see [`smf`](@ref) and
[`invsmf`](@ref)), with `StdUniform` as standard measure. Truncations to
either tail keep their precision (`truncate(StdNormal(), 9..Inf)` has the
mass `1.1e-19`), and so do variates close to the lower end, if `μ` has
accurate interval masses (standard measures and wrapped normal and
logistic distributions do). Transports don't resolve the variates within a
truncation whose uniform variates round to one. See `restrict` for
restrictions to sets given by predicates.
"""
function Base.truncate(μ::AbstractMeasure, a_b::AbstractInterval)
    _check_truncation_parent(μ, _static_ndims(μ))
    _atom_free(rootmeasure(μ)) || _throw_truncation_parent(μ, "a root measure without atoms")
    lo, hi = endpoints(a_b)
    _check_truncation_bounds(dynamic(lo <= hi), a_b)
    return TruncatedMeasure(μ, lo, hi, _checked_smf(μ, smf(μ, lo, hi)))
end

@inline _check_truncation_parent(::Any, ::StaticInteger{0}) = nothing
@inline _check_truncation_parent(μ, ::Any) = _throw_truncation_parent(μ, "scalar variates")
@noinline _throw_truncation_parent(μ, what) =
    throw(ArgumentError("Truncation requires measures with $what, but got a measure of type $(nameof(typeof(μ)))"))

@inline _check_truncation_bounds(ordered::Bool, a_b) = ordered || _throw_truncation_bounds(a_b)
# Traced or device values can't be checked:
@inline _check_truncation_bounds(::Any, a_b) = nothing
@noinline _throw_truncation_bounds(a_b) = throw(ArgumentError("Can't truncate to the empty interval $a_b"))

@inline _checked_smf(μ, F) = F
@noinline _checked_smf(μ, ::NoSMF) =
    throw(ArgumentError("Truncation requires a Stieltjes measure function, measures of type $(nameof(typeof(μ))) have none"))

@inline mspace_elsize(::TruncatedMeasure) = ()
@inline mspace_flatsize(::TruncatedMeasure) = ()
@inline mspace_flatsize(::Type{<:TruncatedMeasure}) = ()
@inline getdof(::TruncatedMeasure) = static(1)
@inline fast_dof(::TruncatedMeasure) = static(1)
@inline preferred_stdmeasure(::Type{<:TruncatedMeasure}) = StdUniform

massof(t::TruncatedMeasure) = t.mass

basemeasure(t::TruncatedMeasure) = basemeasure(t.parent)

@propagate_inbounds checked_arg(t::TruncatedMeasure, x) = checked_arg(t.parent, x)

@inline _in_truncation(t::TruncatedMeasure, x) = (t.lo <= x) & (x <= t.hi)

@inline insupportof(t::TruncatedMeasure, x) = _insupport_and(_in_truncation(t, x), insupportof(t.parent, x))

@inline function logdensityof_impl(t::TruncatedMeasure, x)
    _checksupport(_in_truncation(t, x), _dynamic_logd(logdensityof_impl(t.parent, x), x))
end

@inline logdensity_def(t::TruncatedMeasure, x) = logdensity_def(t.parent, x)

# The Stieltjes measure function of the truncation is the mass of the parent
# above the lower end, its complement the mass below the upper end:
@inline _clamp_to(t::TruncatedMeasure, x) = clamp(x, t.lo, t.hi)
smf(t::TruncatedMeasure, x) = smf(t.parent, t.lo, _clamp_to(t, x))
csmf(t::TruncatedMeasure, x) = smf(t.parent, _clamp_to(t, x), t.hi)
smf(t::TruncatedMeasure, a, x) = smf(t.parent, _clamp_to(t, a), _clamp_to(t, x))

# The inverses are NaN outside of the range `[0, mass]` of the Stieltjes
# measure functions, within it the arguments passed to the parent and the
# results stay within the truncation despite rounding:
invsmf(t::TruncatedMeasure, p) = _truncation_inverse(t, p, invsmf(t.parent, t.lo, _in_smf_range(t, p)))
invcsmf(t::TruncatedMeasure, p) = _truncation_inverse(t, p, invsmf(t.parent, t.hi, -_in_smf_range(t, p)))

@inline _in_smf_range(t::TruncatedMeasure, p) = clamp(p, zero(t.mass), t.mass)

@inline function _truncation_inverse(t::TruncatedMeasure, p, x)
    x_t = clamp(x, t.lo, t.hi)
    ifelse((zero(p) <= p) & (p <= t.mass), x_t, oftype(x_t, NaN))
end
@noinline _truncation_inverse(t::TruncatedMeasure, p, ::NoSMFInverse) =
    throw(ArgumentError("Measures of type $(nameof(typeof(t.parent))) have no inverse Stieltjes measure function"))

# The conditional Stieltjes measure function maps to the unit interval:
@inline transport_to_std(::Type{StdUniform}, t::TruncatedMeasure, x) = _nan_outside(t, x, smf(t, x) / t.mass)
@inline transport_from_std(::Type{StdUniform}, t::TruncatedMeasure, u) =
    _nan_outside(StdUniform(), u, invsmf(t, _unit_interior(u) * t.mass))

Pretty.quoteof(t::TruncatedMeasure) = :(truncate($(Pretty.quoteof(t.parent)), $(t.lo) .. $(t.hi)))
