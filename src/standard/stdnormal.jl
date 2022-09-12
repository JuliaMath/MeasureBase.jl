using SpecialFunctions: erfc, erfcinv
using IrrationalConstants: invsqrt2

struct StdNormal <: StdMeasure end

export StdNormal

@inline insupport(d::StdNormal, x) = true

@inline logdensity_def(::StdNormal, x) = -x^2 / 2
@inline basemeasure(::StdNormal) = WeightedMeasure(static(-0.5 * log2π), LebesgueMeasure())

@inline getdof(::StdNormal) = static(1)

@inline transport_def(::StdUniform, μ::StdNormal, x) = erfc(-z * invsqrt2) / 2
@inline transport_def(::StdNormal, μ::StdUniform, x) = -erfcinv(2*p) * sqrt2

@inline Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdNormal) where {T} = randn(rng, T)


@inline StdMeasure(::typeof(randn)) = StdNormal()
