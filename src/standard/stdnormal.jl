using SpecialFunctions: erfc, erfcinv, logerfc
using IrrationalConstants: invsqrt2, log2π, logtwo, sqrtπ

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
invlogsmf(::StdNormal, logp) = _nan_above_zero(logp, _norminvlogcdf(min(logp, zero(logp))))
invlogcsmf(μ::StdNormal, logp) = -invlogsmf(μ, logp)

# `logerfc` in the far upper tail from the asymptotic series, which is
# exact to the float precision beyond the cut, for implementations of
# `logerfc` that underflow there (Reactant):
@inline function _logerfc(x)
    T = float(real_numtype(typeof(x)))
    cut = sqrt(-log(floatmin(T))) - one(T)
    y = max(x, cut)
    s = inv(2 * y^2)
    series = 1 - s * (1 - 3s * (1 - 5s * (1 - 7s * (1 - 9s * (1 - 11s)))))
    return ifelse(x < cut, logerfc(x), log(series) - y^2 - log(y * T(sqrtπ)))
end

# The standard normal quantile of `exp(lp)`, also where `exp(lp)` underflows:
# from the quantile below the median and the quantile of the complement
# above, in the far tail from the asymptotic `lp ≈ -x^2/2 - log(-x √(2π))`,
# refined by two Newton steps on the log-cdf below the median, branch-free.
# Below `-floatmax/4` the square of the quantile would overflow.
@inline function _norminvlogcdf(lp)
    T = float(real_numtype(typeof(lp)))
    l = max(lp, -floatmax(T) / 4)
    lower = l < -T(logtwo)
    x_q = ifelse(lower, Φinv(_unit_interior(exp(l))), -Φinv(_unit_interior(-expm1(l))))
    t = max(-2 * l, one(T))
    x_tail = -sqrt(max(t - log(t) - T(log2π), zero(T)))
    x = ifelse(l < log(floatmin(T)), x_tail, x_q)
    for _ in 1:2
        x -= (_normlogcdf(x) - l) * exp(_normlogcdf(x) + (x^2 + T(log2π)) / 2)
    end
    return ifelse(lp == -Inf, oftype(x, -Inf), ifelse(lower, x, x_q))
end

smf(::StdNormal) = Φ

transport_def(::StdNormal, μ::StdUniform, p) = _nan_outside(μ, p, Φinv(_unit_interior(p)))
transport_def(::StdUniform, ::StdNormal, x) = Φ(x)
