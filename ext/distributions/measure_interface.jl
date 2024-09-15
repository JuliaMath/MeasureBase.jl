# This file is a part of DistributionMeasures.jl, licensed under the MIT License (MIT).

@inline MeasureBase.logdensity_def(d::Distribution, x) = DensityInterface.logdensityof(d, x)
@inline MeasureBase.unsafe_logdensityof(d::Distribution, x) = DensityInterface.logdensityof(d, x)

@inline MeasureBase.insupport(d::Distribution, x) = Distributions.insupport(d, x)

@inline MeasureBase.basemeasure(d::Distribution{<:ArrayLikeVariate{0},<:Continuous}) = Lebesgue()
@inline MeasureBase.basemeasure(d::Distribution{<:ArrayLikeVariate,<:Continuous}) = Lebesgue()^size(d)
@inline MeasureBase.basemeasure(d::Distribution{<:ArrayLikeVariate{0},<:Discrete}) = Counting()
@inline MeasureBase.basemeasure(d::Distribution{<:ArrayLikeVariate,<:Discrete}) = Counting()^size(d)

@inline MeasureBase.paramnames(d::Distribution) = propertynames(d)
@inline MeasureBase.params(d::Distribution) = NamedTuple{propertynames(d)}(Distributions.params(d))

@inline MeasureBase.testvalue(d::Distribution) = testvalue(basemeasure(d))


@inline MeasureBase.basemeasure(d::Distributions.Poisson) = Counting(MeasureBase.BoundedInts(static(0), static(Inf)))
@inline MeasureBase.basemeasure(d::Distributions.Product{<:Any,<:Distributions.Poisson}) = Counting(MeasureBase.BoundedInts(static(0), static(Inf)))^size(d)


MeasureBase.∫(f, base::Distribution) = MeasureBase.∫(f, convert(AbstractMeasure, base))
