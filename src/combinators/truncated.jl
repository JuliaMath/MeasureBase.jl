export TruncatedMeasure

using LogExpFunctions: logaddexp, logsubexp

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

function _log_interval_mass(μ, a, b)
    m = smf(μ, a, b)
    (a == b || !(m < floatmin(float(typeof(m))))) && return log(m)
    F, G = logsmf(μ, b), logcsmf(μ, a)
    if !(F isa NoSMF) && (G isa NoSMF || F <= G)
        return logsubexp(F, logsmf(μ, a))
    elseif !(G isa NoSMF)
        return logsubexp(G, logcsmf(μ, b))
    end
    log(m)
end
_normalizing_logmass(t::TruncatedMeasure, mass) = _tiny_truncation_mass(t) ? _log_interval_mass(t.parent, t.lo, t.hi) : _logmass(mass)
_sampling_logmass(t::TruncatedMeasure) = _normalizing_logmass(t, t.mass)
logsmf(t::TruncatedMeasure, x) = _log_interval_mass(t.parent, t.lo, _clamp_to(t, x))
logcsmf(t::TruncatedMeasure, x) = _log_interval_mass(t.parent, _clamp_to(t, x), t.hi)

# The inverses are NaN outside of the range `[0, mass]` of the Stieltjes
# measure functions, within it the arguments passed to the parent and the
# results stay within the truncation despite rounding:
invsmf(t::TruncatedMeasure, p) = _truncation_inverse(t, p, invsmf(t.parent, t.lo, _in_smf_range(t, p)))
invcsmf(t::TruncatedMeasure, p) = _truncation_inverse(t, p, invsmf(t.parent, t.hi, -_in_smf_range(t, p)))

function _truncated_invlog(t, lp, upper)
    if !isfinite(lp) || lp >= log(floatmin(float(typeof(lp))))
        return upper ? invcsmf(t, exp(lp)) : invsmf(t, exp(lp))
    end
    lp > _log_interval_mass(t.parent, t.lo, t.hi) && return oftype(float(lp), NaN)
    edge = upper ? t.hi : t.lo
    F, G = logsmf(t.parent, edge), logcsmf(t.parent, edge)
    F isa NoSMF && G isa NoSMF && return upper ? invcsmf(t, exp(lp)) : invsmf(t, exp(lp))
    x = if !(F isa NoSMF) && (G isa NoSMF || F <= G)
        invlogsmf(t.parent, upper ? logsubexp(F, lp) : logaddexp(F, lp))
    else
        invlogcsmf(t.parent, upper ? logaddexp(G, lp) : logsubexp(G, lp))
    end
    clamp(x, t.lo, t.hi)
end
invlogsmf(t::TruncatedMeasure, lp) = _truncated_invlog(t, lp, false)
invlogcsmf(t::TruncatedMeasure, lp) = _truncated_invlog(t, lp, true)

@inline _in_smf_range(t::TruncatedMeasure, p) = clamp(p, zero(t.mass), t.mass)

@inline function _truncation_inverse(t::TruncatedMeasure, p, x)
    x_t = clamp(x, t.lo, t.hi)
    ifelse((zero(p) <= p) & (p <= t.mass), x_t, oftype(x_t, NaN))
end
@noinline _truncation_inverse(t::TruncatedMeasure, p, ::NoSMFInverse) =
    throw(ArgumentError("Measures of type $(nameof(typeof(t.parent))) have no inverse Stieltjes measure function"))

# Scale interval probabilities in log space when the retained mass is subnormal.
_tiny_truncation_mass(t) = t.mass < floatmin(float(typeof(t.mass)))
function _truncated_probability(t, x, upper)
    if _tiny_truncation_mass(t)
        y = _clamp_to(t, x)
        p = upper ? _truncated_interval_probability(t, y, t.hi) : _truncated_interval_probability(t, t.lo, y)
        return oftype(float(y) + t.mass, p)
    end
    return (upper ? csmf(t, x) : smf(t, x)) / t.mass
end
function _truncated_interval_probability(t, a, b)
    F, G = logsmf(t.parent, t.hi), logcsmf(t.parent, t.lo)
    if !(F isa NoSMF) && (G isa NoSMF || F <= G)
        high, low = logsmf(t.parent, b), logsmf(t.parent, a)
    elseif !(G isa NoSMF)
        high, low = logcsmf(t.parent, a), logcsmf(t.parent, b)
    else
        return smf(t.parent, a, b) / t.mass
    end
    p = exp((high - _sampling_logmass(t)) + log(-expm1(min(low - high, zero(high)))))
    return ifelse(high == low, zero(p), p)
end
function _truncated_quantile(t, p, upper)
    if _tiny_truncation_mass(t)
        u = clamp(p, zero(p), one(p))
        left, right = upper ? (u, 1 - u) : (1 - u, u)
        F, G = logsmf(t.parent, t.hi), logcsmf(t.parent, t.lo)
        # Interpolate the endpoint probabilities on the better-resolved tail.
        if !(F isa NoSMF) && (G isa NoSMF || F <= G)
            lp = F + log(right + left * exp(logsmf(t.parent, t.lo) - F))
            x = invlogsmf(t.parent, lp)
        elseif !(G isa NoSMF)
            lp = G + log(left + right * exp(logcsmf(t.parent, t.hi) - G))
            x = invlogcsmf(t.parent, lp)
        else
            return upper ? invcsmf(t, p * t.mass) : invsmf(t, p * t.mass)
        end
        return _nan_outside(StdUniform(), p, _clamp_to(t, x))
    end
    return upper ? invcsmf(t, p * t.mass) : invsmf(t, p * t.mass)
end

smf(n::NormalizedMeasure{<:TruncatedMeasure}, x) = _truncated_probability(n.base, x, false)
csmf(n::NormalizedMeasure{<:TruncatedMeasure}, x) = _truncated_probability(n.base, x, true)
invsmf(n::NormalizedMeasure{<:TruncatedMeasure}, p) = _truncated_quantile(n.base, p, false)
invcsmf(n::NormalizedMeasure{<:TruncatedMeasure}, p) = _truncated_quantile(n.base, p, true)
function smf(n::NormalizedMeasure{<:TruncatedMeasure}, a, x)
    t = n.base
    _tiny_truncation_mass(t) || return smf(t, a, x) / t.mass
    lo, hi = minmax(_clamp_to(t, a), _clamp_to(t, x))
    return copysign(_truncated_interval_probability(t, lo, hi), x - a)
end
function invsmf(n::NormalizedMeasure{<:TruncatedMeasure}, a, p)
    t = n.base
    _tiny_truncation_mass(t) || return invsmf(t, a, p * t.mass)
    return p < 0 ? invcsmf(n, csmf(n, a) - p) : invsmf(n, smf(n, a) + p)
end

# The conditional Stieltjes measure function maps to the unit interval:
@inline transport_to_std(::Type{StdUniform}, t::TruncatedMeasure, x) = _nan_outside(t, x, _truncated_probability(t, x, false))
@inline transport_from_std(::Type{StdUniform}, t::TruncatedMeasure, u) =
    _nan_outside(StdUniform(), u, _truncated_quantile(t, _unit_interior(u), false))

Pretty.quoteof(t::TruncatedMeasure) = :(truncate($(Pretty.quoteof(t.parent)), $(t.lo) .. $(t.hi)))
