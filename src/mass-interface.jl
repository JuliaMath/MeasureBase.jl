import LinearAlgebra: normalize

abstract type MeasureMass end

struct KnownMass{M,MASS} <: MeasureMass
    measure::M
    mass::MASS
end

struct UnknownFiniteMass{M} <: MeasureMass
    measure::M
end

struct UnknownMass{M} <: MeasureMass
    measure::M
end

massof(m::AbstractMeasure) = UnknownMass(m)

struct NormalizedMeasure{P,M} <: AbstractMeasure
    parent::P
    parent_mass::M
end

massof(m::NormalizedMeasure) = KnownMass(m, static(1))

normalize(m::AbstractMeasure) = normalize(m, massof(m))
normalize(m::AbstractMeasure, mass) = NormalizedMeasure(m, mass)

function normalize(m::AbstractMeasure, mass::KnownMass)
    @assert m == mass.measure
    isinf(m.mass) && error("Measure cannot be normalized: $m")
    inv(mass.mass) * m
end

isnormalized(m::AbstractMeasure) = isnormalized(massof(m))
isnormalized(m::KnownMass) = m.mass == 1
isnormalized(::MeasureMass) = false
isnormalized(::NormalizedMeasure) = true
