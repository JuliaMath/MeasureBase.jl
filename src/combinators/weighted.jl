export WeightedMeasure, AbstractWeightedMeasure

"""
    struct WeightedMeasure{R,M} <: AbstractMeasure
        logweight :: R
        base :: M
    end
"""

abstract type AbstractWeightedMeasure <: AbstractMeasure end

# By default the weight for all measure is 1
_logweight(::AbstractMeasure) = 0

@inline function logdensity_def(d::AbstractWeightedMeasure, _)
    d.logweight
end

function Base.rand(rng::AbstractRNG, ::Type{T}, μ::AbstractWeightedMeasure) where {T}
    rand(rng, T, basemeasure(μ))
end

testvalue(::Type{T}, μ::AbstractWeightedMeasure) where {T} = testvalue(T, basemeasure(μ))

###############################################################################

struct WeightedMeasure{R,M} <: AbstractWeightedMeasure
    logweight::R
    base::M
end

massof(w::WeightedMeasure) = exp(w.logweight) * massof(w.base)

_logweight(μ::WeightedMeasure) = μ.logweight
basemeasure(μ::AbstractWeightedMeasure) = μ.base

function Pretty.tile(d::WeightedMeasure)
    weight = round(dynamic(exp(d.logweight)), sigdigits = 4)
    Pretty.pair_layout(Pretty.tile(weight), Pretty.tile(d.base), sep = " * ")
end

function Base.:*(k::T, m::AbstractMeasure) where {T<:Number}
    logk = log(k)
    return weightedmeasure(logk, m)
end

Base.:*(m::AbstractMeasure, k::Real) = k * m

gentype(μ::WeightedMeasure) = gentype(μ.base)

insupport(μ::WeightedMeasure, x) = insupport(μ.base, x)

# TODO: Transports must preserve mass
transport_origin(ν::WeightedMeasure) = ν.base

to_origin(w::WeightedMeasure, y) = y
from_origin(w::WeightedMeasure, x) = x
