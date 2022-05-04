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

Base.length(μ::AbstractWeightedMeasure) = length(basemeasure(μ))

@inline function logdensity_def(d::AbstractWeightedMeasure, _)
    d.logweight
end

function Base.rand(rng::AbstractRNG, ::Type{T}, μ::AbstractWeightedMeasure) where {T}
    rand(rng, T, basemeasure(μ))
end

###############################################################################

struct WeightedMeasure{R,M} <: AbstractWeightedMeasure
    logweight::R
    base::M
end

_logweight(μ::WeightedMeasure) = μ.logweight
basemeasure(μ::AbstractWeightedMeasure) = μ.base

function Base.show(io::IO, μ::WeightedMeasure)
    io = IOContext(io, :compact => true)
    print(io, exp(μ.logweight), " * ", μ.base)
end

function Base.show_unquoted(io::IO, μ::WeightedMeasure, indent::Int, prec::Int)
    io = IOContext(io, :compact => true)
    if Base.operator_precedence(:*) ≤ prec
        print(io, "(")
        show(io, μ)
        print(io, ")")
    else
        show(io, μ)
    end
    return nothing
end

function Base.:*(k::T, m::AbstractMeasure) where {T<:Number}
    logk = log(k)
    return weightedmeasure(logk, m)
end

Base.:*(m::AbstractMeasure, k::Real) = k * m

≪(::M, ::WeightedMeasure{R,M}) where {R,M} = true
≪(::WeightedMeasure{R,M}, ::M) where {R,M} = true

gentype(μ::WeightedMeasure) = gentype(μ.base)

insupport(μ::WeightedMeasure, x) = insupport(μ.base, x)