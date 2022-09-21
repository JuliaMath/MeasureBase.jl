import LinearAlgebra: normalize

import Base

abstract type AbstractUnknownMass <: Number end

struct UnknownFiniteMass <: AbstractUnknownMass
end

struct UnknownMass <: AbstractUnknownMass
end

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

massof(m::AbstractMeasure) = UnknownMass(m)

logmassof(m::AbstractMeasure) = log(massof(m))

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
isnormalized(x, p::Real=2) = isone(norm(x, p))

isone(::AbstractUnknownMass) = false