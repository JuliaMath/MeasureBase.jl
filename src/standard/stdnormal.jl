using SpecialFunctions

struct StdNormal <: StdMeasure end

export StdNormal

@inline insupport(d::StdNormal, x) = true

@inline logdensity_def(::StdNormal, x) = -x^2 / 2
@inline basemeasure(::StdNormal) = WeightedMeasure(static(-0.5 * log2π), LebesgueMeasure())

@inline getdof(::StdNormal) = static(1)

@inline transport_def(::StdUniform, μ::StdNormal, x) = StatsFuns.normcdf(x)
@inline transport_def(::StdNormal, μ::StdUniform, x) = StatsFuns.norminvcdf(x)

@inline Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdNormal) where {T} = randn(rng, T)


@inline StdMeasure(::typeof(randn)) = StdNormal()