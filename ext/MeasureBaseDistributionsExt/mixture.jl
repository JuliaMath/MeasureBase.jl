# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Mixtures become superpositions, except continuous univariate mixtures,
# which keep the cdf/quantile-based transport of wrapped distributions:
function MeasureBase.AbstractMeasure(d::Distributions.AbstractMixtureModel)
    superpose(map((w, c) -> w * asmeasure(c), Distributions.probs(d), Distributions.components(d)))
end
MeasureBase.AbstractMeasure(d::Distributions.AbstractMixtureModel{Univariate,Continuous}) = AsMeasure{typeof(d)}(d)


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
