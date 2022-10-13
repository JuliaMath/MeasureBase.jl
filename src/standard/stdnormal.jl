using SpecialFunctions: erfc, erfcinv
using IrrationalConstants: invsqrt2

struct StdNormal <: StdMeasure end

export StdNormal

@inline insupport(d::StdNormal, x) = true

@inline logdensity_def(::StdNormal, x) = -x^2 / 2
@inline basemeasure(::StdNormal) = WeightedMeasure(static(-0.5 * log2π), LebesgueBase())

@inline getdof(::StdNormal) = static(1)

@inline Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdNormal) where {T} = randn(rng, T)

Φ(z) = erfc(-z * invsqrt2) / 2
Φinv(p) = -erfcinv(2 * p) * sqrt2

smf(::StdNormal, x) = Φ(x)
invsmf(::StdNormal, p) = Φinv(p)

transport_def(::StdNormal, ::StdUniform, p) = invsmf(StdNormal(), p)
transport_def(::StdUniform, ::StdNormal, x) = smf(StdNormal(), x)
