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
`dμ/dν`, also called the _Radom-Nikodym derivative_:
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
`log(dμ/dν)`, also called the _Radom-Nikodym derivative_:
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

A `DensityMeasure` is a measure defined by a density or log-density with
respect to some other "base" measure.

Users should not instantiate `DensityMeasure` directly, but should instead
call `mintegral(f, base)` (if `f` is a density function or
`DensityInterface.IsDensity` object) or `mintegral_exp(f, base)` (if `f`
is a log-density function).
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
    result = Pretty.literal("mintegrate(")
    result *= Pretty.pair_layout(Pretty.tile(μ.f), Pretty.tile(μ.base); sep = ", ")
    result *= Pretty.literal(")")
end

basemeasure(μ::DensityMeasure) = μ.base

logdensity_def(μ::DensityMeasure, x) = logdensityof(μ.f, x)

density_def(μ::DensityMeasure, x) = densityof(μ.f, x)

localmeasure(μ::DensityMeasure, x) = DensityMeasure(μ.f, localmeasure(μ.base, x))

@doc raw"""
    mintegrate(f, μ::AbstractMeasure)::AbstractMeasure

Returns a new measure that represents the indefinite
[integral](https://en.wikipedia.org/wiki/Radon%E2%80%93Nikodym_theorem)
of `f` with respect to `μ`.

`ν = mintegrate(f, μ)` generates a measure `ν` that has the mathematical
interpretation

math```
\nu(A) = \int_A f(a) \, \rm{d}\mu(a)
```
"""
function mintegrate end
export mintegrate

mintegrate(f, μ::AbstractMeasure) = _mintegrate_impl(f, μ, DensityKind(f))

_mintegrate_impl(f, μ, ::IsDensity) = DensityMeasure(f, μ)
function _mintegrate_impl(f, μ, ::HasDensity)
    throw(
        ArgumentError(
            "`mintegrate(f, mu)` requires `DensityKind(f)` to be `IsDensity()` or `NoDensity()`.",
        ),
    )
end
_mintegrate_impl(f, μ, ::NoDensity) = DensityMeasure(funcdensity(f), μ)

@doc raw"""
    mintegrate_exp(log_f, μ::AbstractMeasure)

Given a function `log_f` that semantically represents the log of a function
`f`, `mintegrate` returns a new measure that represents the indefinite
[integral](https://en.wikipedia.org/wiki/Radon%E2%80%93Nikodym_theorem)
of `f` with respect to `μ`.

`ν = mintegrate_exp(log_f, μ)` generates a measure `ν` that has the
mathematical interpretation

math```
\nu(A) = \int_A e^{log(f(a))} \, \rm{d}\mu(a) = \int_A f(a) \, \rm{d}\mu(a)
```

Note that `exp(log_f(...))` is usually not run explicitly, calculations that
involve the resulting measure are typically performed in log-space,
internally.
"""
function mintegrate_exp end
export mintegrate_exp

function mintegrate_exp(log_f, μ::AbstractMeasure)
    _mintegrate_exp_impl(log_f, μ, DensityKind(log_f))
end

function _mintegrate_exp_impl(log_f, μ, ::IsDensity)
    throw(
        ArgumentError(
            "`mintegrate_exp(log_f, μ)` is not valid when `DensityKind(log_f) == IsDensity()`. Use `mintegrate(log_f, μ)` instead.",
        ),
    )
end
function _mintegrate_exp_impl(log_f, μ, ::HasDensity)
    throw(
        ArgumentError(
            "`mintegrate_exp(log_f, μ)` is not valid when `DensityKind(log_f) == HasDensity()`.",
        ),
    )
end
_mintegrate_exp_impl(log_f, μ, ::NoDensity) = DensityMeasure(logfuncdensity(log_f), μ)
