export ↑

struct PowerWeightedMeasure{M,A} <: AbstractMeasure
    parent::M
    exponent::A
end

logdensity_def(d::PowerWeightedMeasure, x) = d.exponent * logdensity_def(d.parent, x)

basemeasure(d::PowerWeightedMeasure, x) = basemeasure(d.parent, x)

basemeasure(d::PowerWeightedMeasure) = basemeasure(d.parent) ↑ d.exponent

function powerweightedmeasure(d, α)
    isone(α) && return d
    PowerWeightedMeasure(d, α)
end

(d::AbstractMeasure) ↑ α = powerweightedmeasure(d, α)

insupport(d::PowerWeightedMeasure, x) = insupport(d.parent, x)
