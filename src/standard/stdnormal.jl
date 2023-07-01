using SpecialFunctions: erfc, erfcinv
using IrrationalConstants: invsqrt2

"""
    StdNormal <: StdMeasure

Represents the standard (mean of zero, variance of one)
[normal](https://en.wikipedia.org/wiki/Normal_distribution) probability measure.

See [`StdMeasure`](@ref) for the semantics of standard measures in the
context of MeasureBase.
"""
struct StdNormal <: StdMeasure end
export StdNormal

@inline insupport(::StdNormal, x) = true

@inline logdensity_def(::StdNormal, x) = -x^2 / 2
@inline basemeasure(::StdNormal) = WeightedMeasure(static(-0.5 * log2π), LebesgueBase())

@inline getdof(::StdNormal) = static(1)

@inline Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdNormal) where {T} = randn(rng, T)

Φ(z) = erfc(-z * invsqrt2) / 2
Φinv(p) = -erfcinv(2 * p) * sqrt2

InverseFunctions.inverse(::typeof(Φ)) = Φinv
InverseFunctions.inverse(::typeof(Φinv)) = Φ

smf(::StdNormal, x) = Φ(x)
invsmf(::StdNormal, p) = Φinv(p)

smf(::StdNormal) = Φ
invsmf(::StdNormal) = Φinv

transport_def(::StdNormal, ::StdUniform, p) = Φinv(p)
transport_def(::StdUniform, ::StdNormal, x) = Φ(x)
