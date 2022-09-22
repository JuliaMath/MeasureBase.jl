export ↑

struct PowerWeightedMeasure{M,A} <: AbstractMeasure
    parent::M
    exponent::A
end

logdensity_def(d::PowerWeightedMeasure, x) = d.exponent * logdensity_def(d.parent, x)

basemeasure(d::PowerWeightedMeasure, x) = basemeasure(d.parent, x)↑d.exponent

basemeasure(d::PowerWeightedMeasure) = basemeasure(d.parent)↑d.exponent

function powerweightedmeasure(d, α)
    isone(α) && return d
    PowerWeightedMeasure(d, α)
end

(d::AbstractMeasure)↑α = powerweightedmeasure(d, α)

insupport(d::PowerWeightedMeasure, x) = insupport(d.parent, x)

function Base.show(io::IO, d::PowerWeightedMeasure)
    print(io, d.parent, " ↑ ", d.exponent)
end

function powerweightedmeasure(d::PowerWeightedMeasure, α)
    powerweightedmeasure(d.parent, α * d.exponent)
end

function powerweightedmeasure(d::WeightedMeasure, α)
    weightedmeasure(α * d.logweight, powerweightedmeasure(d.base, α))
end

function Pretty.tile(d::PowerWeightedMeasure)
    Pretty.pair_layout(Pretty.tile(d.parent), Pretty.tile(d.exponent), sep = " ↑ ")
end

massof(m::PowerWeightedMeasure) = massof(m.parent)^m.exponent
