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

Instead of calling this directly, users should call `logdensity_rel(μ, ν)`.
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
call `mintegrate(f, base)` (if `f` is a density function or
`DensityInterface.IsDensity` object) or `mintegrate_exp(f, base)` (if `f`
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


"""
    MeasureBase.as_integrand(f)
    MeasureBase.as_integrand(density)

Make `f` or `density` (more) suitable as an integrand for
[`mintegrate`](@ref).

`mintegrate(obj, μ::AbstractMeasure)` automatically calls
`as_integrand(obj)` internally.

If a density is passed, it must implement the DensityInterface API.

By default just returns `f` resp. `density`, but may be specialized for
functions and densities that can profit from conversion to a form optimized
for use in `mintegrate`.

See also [`MeasureBase.as_likelihood`](@ref).
"""
function as_integrand end

@inline as_integrand(obj) = _as_integrand_default_impl(obj, DensityKind(obj))

@inline _as_integrand_default_impl(f, ::NoDensity) = funcdensity(f)

@inline _as_integrand_default_impl(density, ::IsDensity) = density

function _as_integrand_default_impl(obj, ::HasDensity)
    throw(
        ArgumentError(
            "`MeasureBase.as_integrand(obj)` requires `DensityKind(obj)` to be `IsDensity()` or `NoDensity()`.",
        ),
    )
end


@doc raw"""
    mintegrate(f, μ::AbstractMeasure)::AbstractMeasure
    mintegrate(density, μ::AbstractMeasure)::AbstractMeasure

Returns a new measure that represents the indefinite
[integral](https://en.wikipedia.org/wiki/Radon%E2%80%93Nikodym_theorem)
of `f` with respect to `μ`.

If a density is passed, it must implement the DensityInterface API.

`ν = mintegrate(f, μ)` generates a measure `ν` that has the mathematical
interpretation

```math
\nu(A) = \int_A f(a) \, \rm{d}\mu(a)
```
"""
function mintegrate end
export mintegrate

@inline mintegrate(obj, μ::AbstractMeasure) = DensityMeasure(as_integrand(obj), μ)


"""
    MeasureBase.as_integrand_exp(log_f)

Convert the logarithm of an integrand to an integrand.

See also [`MeasureBase.as_integrand`](@ref).
"""
function as_integrand_exp end

@inline as_integrand_exp(log_f) = _as_integrand_exp_default_impl(log_f, DensityKind(log_f))

@inline _as_integrand_exp_default_impl(log_f, ::NoDensity) = logfuncdensity(log_f)

function _as_integrand_exp_default_impl(log_f, ::IsDensity)
    throw(
        ArgumentError(
            "`as_integrand_exp(log_f)` is not valid when `DensityKind(log_f) == IsDensity()`. Use `as_integrand(log_f)` instead.",
        ),
    )
end

function _as_integrand_exp_default_impl(log_f, ::HasDensity)
    throw(
        ArgumentError(
            "`as_integrand_exp(log_f)` is not valid when `DensityKind(log_f) == HasDensity()`.",
        ),
    )
end


@doc raw"""
    mintegrate_exp(log_f, μ::AbstractMeasure)

Given a function `log_f` that semantically represents the log of a function
`f`, `mintegrate_exp` returns a new measure that represents the indefinite
[integral](https://en.wikipedia.org/wiki/Radon%E2%80%93Nikodym_theorem)
of `f` with respect to `μ`.

`ν = mintegrate_exp(log_f, μ)` generates a measure `ν` that has the
mathematical interpretation

```math
\nu(A) = \int_A e^{log(f(a))} \, \rm{d}\mu(a) = \int_A f(a) \, \rm{d}\mu(a)
```

Note that `exp(log_f(...))` is usually not run explicitly, calculations that
involve the resulting measure are typically performed in log-space,
internally.
"""
function mintegrate_exp end
export mintegrate_exp

mintegrate_exp(log_f, μ::AbstractMeasure) = DensityMeasure(as_integrand_exp(log_f), μ)


basemeasure(μ::DensityMeasure) = μ.base

logdensity_def(μ::DensityMeasure, x) = logdensityof(μ.f, x)

density_def(μ::DensityMeasure, x) = densityof(μ.f, x)

function logdensityof_impl(μ::DensityMeasure, x::Any)
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
