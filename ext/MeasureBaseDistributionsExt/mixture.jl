# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

function MeasureBase.AbstractMeasure(d::Distributions.AbstractMixtureModel)
    superpose(map((w, c) -> w * asmeasure(c), Distributions.probs(d), Distributions.components(d)))
end

function AsMeasure{D}(::D) where {D<:Distributions.AbstractMixtureModel}
    throw(ArgumentError("Don't wrap Distributions.AbstractMixtureModel into MeasureBase.AsMeasure, use asmeasure to convert instead."))
end


const _MixtureMeasure = SuperpositionMeasure{
    <:Union{Tuple{Vararg{WeightedMeasure}},AbstractVector{<:WeightedMeasure}},
}

_mixture_component(m::AsMeasure{<:Distribution}) = m.obj

function Distributions.Distribution(m::_MixtureMeasure)
    components = map(c -> _mixture_component(c.base), collect(values(m.components)))
    prior = map(c -> exp(c.logweight), collect(values(m.components)))
    Distributions.MixtureModel(components, prior)
end

Base.convert(::Type{Distribution}, m::_MixtureMeasure) = Distributions.Distribution(m)
