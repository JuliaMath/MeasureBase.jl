struct StdNormal <: MeasureBase.StdMeasure end

export StdNormal

@inline MeasureBase.insupport(d::StdNormal, x) = true

@inline MeasureBase.logdensity_def(::StdNormal, x) = -x^2 / 2
@inline MeasureBase.basemeasure(::StdNormal) = WeightedMeasure(static(-0.5 * log2π), Lebesgue(ℝ))

@inline MeasureBase.getdof(::StdNormal) = static(1)

@inline MeasureBase.transport_def(::StdUniform, μ::StdNormal, x) = StatsFuns.normcdf(x)
@inline MeasureBase.transport_def(::StdNormal, μ::StdUniform, x) = StatsFuns.norminvcdf(x)

@inline Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdNormal) where {T} = randn(rng, T)


@inline MeasureBase.StdMeasure(::typeof(randn)) = StdNormal()