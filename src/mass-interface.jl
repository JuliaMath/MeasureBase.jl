import LinearAlgebra: normalize

import Base

"""
    abstract type MeasureBase.AbstractUnknownMass <: Number

Supertype of [`UnknownFiniteMass`](@ref) and [`UnknownMass`](@ref), the
results of [`massof`](@ref) when the mass of a measure is not known.
"""
abstract type AbstractUnknownMass <: Number end

"""
    struct MeasureBase.UnknownFiniteMass <: AbstractUnknownMass

The mass of a measure is finite but unknown. See [`massof`](@ref).
"""
struct UnknownFiniteMass <: AbstractUnknownMass end

"""
    struct MeasureBase.UnknownMass <: AbstractUnknownMass

The mass of a measure is unknown, possibly infinite. See [`massof`](@ref).
"""
struct UnknownMass <: AbstractUnknownMass end

for T in (:UnknownFiniteMass, :UnknownMass)
    @eval begin
        Base.:+(::$T, ::$T) = $T()
        Base.:+(x::Real, ::$T) = isfinite(dynamic(x)) ? $T() : UnknownMass()
        Base.:+(::$T, x::Real) = isfinite(dynamic(x)) ? $T() : UnknownMass()
        Base.:*(::$T, ::$T) = $T()
        Base.:*(x::Real, ::$T) = isfinite(dynamic(x)) ? $T() : UnknownMass()
        Base.:*(::$T, x::Real) = isfinite(dynamic(x)) ? $T() : UnknownMass()
        Base.:^(::$T, k::Real) = isfinite(k) ? $T() : UnknownMass()
        # Disambiguation:
        Base.:^(::$T, k::Integer) = isfinite(k) ? $T() : UnknownMass()
        Base.:^(::$T, k::Rational) = isfinite(k) ? $T() : UnknownMass()
    end
end

for op in (:+, :*)
    let
        U = :UnknownMass
        UF = :UnknownFiniteMass
        @eval begin
            Base.$op(::$U, ::$UF) = $U()
            Base.$op(::$UF, ::$U) = $U()
        end
    end
end

export massof

"""
    massof(m)

Get the _mass_ of a measure - that is, integrate the measure over its support.

----------

    massof(m, dom)

Integrate the measure `m` over the "domain" `dom`. Note that domains are not
defined universally, but may be specific to a given measure. If `m` is
`<:AbstractMeasure`, users can also write `m(dom)`. For new measures, users
should *not* add new "call" methods, but instead extend `MeasureBase.massof`.


For example, for many univariate measures `m` with `rootmeasure(m) ==
LebesgueBase()`, users can call `massof(m, a_b)` where
`a_b::IntervalSets.Interval`.

`massof` often returns a `Real`. But in many cases we may only know the mass is
finite, or we may know nothing at all about it. For these cases, it will return
[`UnknownFiniteMass()`](@ref UnknownFiniteMass) or [`UnknownMass()`](@ref UnknownMass),
respectively. When no `massof` method exists, it defaults to `UnknownMass()`.

`massof` must be type stable. Types of normalized measures should return
`static(1.0)`, [`normalize`](@ref) then returns them as they are.
"""
massof(::AbstractMeasure) = UnknownMass()

# The log of a known mass, densities that involve unknown masses are unknown:
@inline _logmass(mass::Real) = log(dynamic(mass))
@noinline _logmass(::AbstractUnknownMass) =
    throw(ArgumentError("Density depends on a measure of unknown mass"))

export isnormalized

"""
    isnormalized(m::AbstractMeasure)

Checks whether the measure m is normalized, that is, whether `massof(m) == 1`. 

For convenience, we also provide a method on non-measures that only depends on
`norm`.
"""
isnormalized(m::AbstractMeasure) = isone(massof(m))

"""
    isnormalized(x, p::Real=2)

Check whether `norm(x, p) == 1`.
"""
isnormalized(x, p::Real = 2) = isone(norm(x, p))

Base.isone(::AbstractUnknownMass) = false

function massof(m, s)
    _default_massof_impl(m, s, rootmeasure(m))
end

"""
    (m::AbstractMeasure)(s)

Convenience method for `massof(m, s)`. To make a user-defined measure callable
in this way, users should add the corresponding `massof` method.
"""
(m::AbstractMeasure)(s) = massof(m, s)

function massof(μ, a_b::AbstractInterval)
    _smf_interval_massof(μ, smf(μ, rightendpoint(a_b)), smf(μ, leftendpoint(a_b)))
end

_smf_interval_massof(μ, smf_r, smf_l) = smf_r - smf_l

function _smf_interval_massof(μ, ::NoSMF, ::NoSMF)
    throw(
        ArgumentError(
            "Can't compute the mass over an interval for a measure of type $(nameof(typeof(μ))), no statistical measure function available",
        ),
    )
end
