# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Bare distributions behave like their measures (some constructors keep
# them, e.g. products of distributions). Their densities are Distributions'
# ones, relative to the base measures of their measures:
@inline MeasureBase.logdensity_def(d::Distribution, x) = DensityInterface.logdensityof(d, x)
@inline MeasureBase.unsafe_logdensityof(d::Distribution, x) = DensityInterface.logdensityof(d, x)
@inline MeasureBase.insupportof(d::Distribution, x) = Distributions.insupport(d, x) & _finite_variate(d, x)
@inline MeasureBase.basemeasure(d::Distribution) = basemeasure(asmeasure(d))
@inline MeasureBase.rootmeasure(d::Distribution) = rootmeasure(asmeasure(d))
@inline MeasureBase.testvalue(d::Distribution) = testvalue(asmeasure(d))
@inline MeasureBase.testvalue(::Type{T}, d::Distribution) where {T} = testvalue(T, asmeasure(d))

MeasureBase.mintegrate(f, base::Distribution) =
    MeasureBase.mintegrate(f, convert(AbstractMeasure, base))


@inline MeasureBase.paramnames(d::Distribution) = propertynames(d)
@inline MeasureBase.paramsof(d::Distribution) = NamedTuple{propertynames(d)}(Distributions.params(d))
