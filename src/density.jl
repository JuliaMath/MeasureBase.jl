###################################################################
# Abstract types and methods

abstract type AbstractDensity <: Function end

@inline DensityKind(::AbstractDensity) = IsDensity()

import DensityInterface

####################################################################################
# Density

"""
    struct Density{M,B} <: AbstractDensity
        μ::M
        base::B
    end

For measures `μ` and `ν`, `Density(μ,ν)` represents the _density function_
`dμ/dν`, also called the _Radon-Nikodym derivative_:
https://en.wikipedia.org/wiki/Radon%E2%80%93Nikodym_theorem#Radon%E2%80%93Nikodym_derivative

Instead of calling this directly, users should call `density_rel(μ, ν)`.
"""
struct Density{M,B} <: AbstractDensity
    μ::M
    base::B
end

Base.:∘(::typeof(log), d::Density) = logdensity_rel(d.μ, d.base)

Base.log(d::Density) = log ∘ d

density_rel(μ, base) = Density(μ, base)

(f::Density)(x) = density_rel(f.μ, f.base, x)

DensityInterface.logfuncdensity(d::Density) = throw(MethodError(logfuncdensity, (d,)))

####################################################################################
# LogDensity

"""
    struct LogDensity{M,B} <: AbstractDensity
        μ::M
        base::B
    end

For measures `μ` and `ν`, `LogDensity(μ,ν)` represents the _log-density function_
`log(dμ/dν)`, also called the _Radon-Nikodym derivative_:
https://en.wikipedia.org/wiki/Radon%E2%80%93Nikodym_theorem#Radon%E2%80%93Nikodym_derivative

Instead of calling this directly, users should call `logdensity_rel(μ, ν)` or
its abbreviated form, `log𝒹(μ,ν)`.
"""
struct LogDensity{M,B} <: AbstractDensity
    μ::M
    base::B
end

Base.:∘(::typeof(exp), d::LogDensity) = density_rel(d.μ, d.base)

Base.exp(d::LogDensity) = exp ∘ d

logdensity_rel(μ, base) = LogDensity(μ, base)

(f::LogDensity)(x) = logdensity_rel(f.μ, f.base, x)

DensityInterface.funcdensity(d::LogDensity) = throw(MethodError(funcdensity, (d,)))

#######################################################################################
# DensityMeasure

"""
    struct DensityMeasure{F,B} <: AbstractDensityMeasure
        density :: F
        base    :: B
    end

A `DensityMeasure` is a measure defined by a density or log-density with respect
to some other "base" measure.

Users should not call `DensityMeasure` directly, but should instead call `∫(f,
base)` (if `f` is a density function or `DensityInterface.IsDensity` object) or
`∫exp(f, base)` (if `f` is a log-density function).
"""
struct DensityMeasure{F,B} <: AbstractMeasure
    f::F
    base::B

    function DensityMeasure(f::F, base::B) where {F,B}
        @assert DensityKind(f) isa IsDensity
        new{F,B}(f, base)
    end
end

@inline function insupport(d::DensityMeasure, x)
    insupport(d.base, x) == true && isfinite(logdensityof(getfield(d, :f), x))
end

function Pretty.tile(μ::DensityMeasure{F,B}) where {F,B}
    result = Pretty.literal("DensityMeasure ∫(")
    result *= Pretty.pair_layout(Pretty.tile(μ.f), Pretty.tile(μ.base); sep = ", ")
    result *= Pretty.literal(")")
end

export ∫

"""
    ∫(f, base::AbstractMeasure)

Define a new measure in terms of a density `f` over some measure `base`.
"""
∫(f, base) = _densitymeasure(f, base, DensityKind(f))

_densitymeasure(f, base, ::IsDensity) = DensityMeasure(f, base)
function _densitymeasure(f, base, ::HasDensity)
    @error "`∫(f, base)` requires `DensityKind(f)` to be `IsDensity()` or `NoDensity()`."
end
_densitymeasure(f, base, ::NoDensity) = DensityMeasure(funcdensity(f), base)

export ∫exp

"""
    ∫exp(f, base::AbstractMeasure)

Define a new measure in terms of a log-density `f` over some measure `base`.
"""
∫exp(f, base) = _logdensitymeasure(f, base, DensityKind(f))

function _logdensitymeasure(f, base, ::IsDensity)
    @error "`∫exp(f, base)` is not valid when `DensityKind(f) == IsDensity()`. Use `∫(f, base)` instead."
end
function _logdensitymeasure(f, base, ::HasDensity)
    @error "`∫exp(f, base)` is not valid when `DensityKind(f) == HasDensity()`."
end
_logdensitymeasure(f, base, ::NoDensity) = DensityMeasure(logfuncdensity(f), base)

basemeasure(μ::DensityMeasure) = μ.base

logdensity_def(μ::DensityMeasure, x) = logdensityof(μ.f, x)

density_def(μ::DensityMeasure, x) = densityof(μ.f, x)

function logdensityof(μ::DensityMeasure, x::Any)
    integrand, μ_base = μ.f, μ.base

    base_logval = logdensityof(μ_base, x)

    T = typeof(base_logval)
    U = logdensityof_rt(integrand, x)
    R = promote_type(T, U)

    # Don't evaluate base measure if integrand is zero or NaN
    if isneginf(base_logval)
        R(-Inf)
    else
        integrand_logval = logdensityof(integrand, x)
        convert(R, integrand_logval + base_logval)::R
    end
end
