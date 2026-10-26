using SpecialFunctions: erfc, erfcinv, logerfc
using IrrationalConstants: invsqrt2, log2π, logtwo

"""
    StdNormal <: StdMeasure

The standard normal measure, the normal distribution with zero mean and
unit variance as a measure.
"""
struct StdNormal <: StdMeasure end

export StdNormal

@inline insupportof(::StdNormal, x) = true

# Integer variates are converted before squaring, which could overflow:
@inline logdensityof_impl(::StdNormal, x) = (-float(x)^2 - log2π) / 2

@inline logdensity_def(::StdNormal, x) = -float(x)^2 / 2
@inline basemeasure(::StdNormal) = WeightedMeasure(static(-0.5 * log2π), LebesgueBase())

@inline rand_impl(ctx::GenContext, ::StdNormal) = randn(ctx)
@inline batched_rand_impl(ctx::GenContext, ::StdNormal, sz::SizeLike) = _randn_bulk(ctx, sz)

Φ(z) = erfc(-z * invsqrt2) / 2
Φinv(p) = -erfcinv(2 * p) * sqrt2

InverseFunctions.inverse(::typeof(Φ)) = Φinv
InverseFunctions.inverse(::typeof(Φinv)) = Φ

smf(::StdNormal, x) = Φ(x)
invsmf(μ::StdNormal, p) = transport_def(μ, StdUniform(), p)
csmf(::StdNormal, x) = _normccdf(x)
invcsmf(μ::StdNormal, p) = -invsmf(μ, p)
# The log-cdf rounds to zero above the median, the log of the complement is
# accurate there:
logsmf(::StdNormal, x) = ifelse(x < zero(x), _normlogcdf(x), log1p(-_normccdf(x)))
logcsmf(::StdNormal, x) = ifelse(x > zero(x), _normlogccdf(x), log1p(-Φ(x)))

smf(::StdNormal) = Φ

transport_def(::StdNormal, μ::StdUniform, p) = _nan_outside(μ, p, Φinv(_unit_interior(p)))
transport_def(::StdUniform, ::StdNormal, x) = Φ(x)
