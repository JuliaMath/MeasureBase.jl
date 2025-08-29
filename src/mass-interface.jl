
import LinearAlgebra: normalize

import Base

abstract type AbstractUnknownMass <: Number end

"""
    struct UnknownFiniteMass <: AbstractUnknownMass end

See `massof`
"""
struct UnknownFiniteMass <: AbstractUnknownMass end

"""
    struct UnknownMass <: AbstractUnknownMass end

See `massof`
"""
struct UnknownMass <: AbstractUnknownMass end

for T in (:UnknownFiniteMass, :UnknownMass)
    @eval begin
        Base.:+(::$T, ::$T) = $T()
        Base.:*(::$T, ::$T) = $T()
        Base.:^(::$T, k::Number) = isfinite(k) ? $T() : UnknownMass()
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

`massof` 

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
`UnknownFiniteMass` or `UnknownMass`, respectively. When no `massof` method
exists, it defaults to `UnknownMass`.
"""
massof(m::AbstractMeasure) = UnknownMass(m)

struct NormalizedMeasure{P,M} <: AbstractMeasure
    parent::P
    parent_mass::M
end

massof(m::NormalizedMeasure) = static(1.0)

normalize(m::AbstractMeasure) = _normalize(m, massof(m))

_normalize(m::AbstractMeasure, mass::AbstractUnknownMass) = NormalizedMeasure(m, mass)

function _normalize(m::AbstractMeasure, mass)
    isinf(mass) && error("Measure cannot be normalized: $m")
    inv(mass) * m
end

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

isone(::AbstractUnknownMass) = false

function massof(m, s)
    _massof(m, s, rootmeasure(m))
end

"""
    (m::AbstractMeasure)(s)

Convenience method for `massof(m, s)`. To make a user-defined measure callable
in this way, users should add the corresponding `massof` method.
"""
(m::AbstractMeasure)(s) = massof(m, s)

massof(μ, a_b::AbstractInterval) = smf(μ, rightendpoint(a_b)) - smf(μ, leftendpoint(a_b))
