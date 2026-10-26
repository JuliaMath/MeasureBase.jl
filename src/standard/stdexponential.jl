"""
    StdExponential <: StdMeasure

The standard exponential measure, the exponential distribution with unit
scale as a measure.
"""
struct StdExponential <: StdMeasure end

export StdExponential

insupportof(::StdExponential, x) = x ≥ zero(x)

@inline function logdensityof_impl(d::StdExponential, x)
    R = float(typeof(x))
    _checksupport(insupportof(d, x), convert(R, -x))
end

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = LebesgueBase()

@inline transport_def(::StdUniform, μ::StdExponential, x) = _nan_outside(μ, x, -expm1(-x))
@inline transport_def(::StdExponential, μ::StdUniform, x) = _nan_outside(μ, x, -log1p(-_unit_interior(x)))

@inline rand_impl(ctx::GenContext, ::StdExponential) = randexp(ctx)
@inline batched_rand_impl(ctx::GenContext, ::StdExponential, sz::SizeLike) = _randexp_bulk(ctx, sz)

smf(::StdExponential, x) = -expm1(-max(x, zero(x)))
invsmf(μ::StdExponential, p) = transport_def(μ, StdUniform(), p)
csmf(::StdExponential, x) = exp(-max(x, zero(x)))
invcsmf(::StdExponential, p) = _nan_outside(StdUniform(), p, -log(_unit_interior(p)))
logsmf(::StdExponential, x) = log1mexp(-max(float(x), zero(float(x))))
logcsmf(::StdExponential, x) = -max(float(x), zero(float(x)))
invlogsmf(::StdExponential, logp) = _nan_above_zero(logp, -log1mexp(min(logp, zero(logp))))
invlogcsmf(::StdExponential, logp) = _nan_above_zero(logp, -float(logp))

# Log-probabilities above zero are invalid arguments:
@inline _nan_above_zero(logp, y) = ifelse(logp > zero(logp), oftype(y, NaN), y)
