"""
    StdLogistic <: StdMeasure

The standard logistic measure, the logistic distribution with zero
location and unit scale as a measure.
"""
struct StdLogistic <: StdMeasure end

export StdLogistic

@inline insupportof(d::StdLogistic, x) = true

@inline logdensityof_impl(::StdLogistic, x) = (u = -abs(x); u - 2 * log1pexp(u))

@inline logdensity_def(::StdLogistic, x) = logdensityof(StdLogistic(), x)
@inline basemeasure(::StdLogistic) = LebesgueBase()

@inline transport_def(::StdUniform, μ::StdLogistic, x) = logistic(x)
@inline transport_def(::StdLogistic, μ::StdUniform, p) = _nan_outside(μ, p, logit(_unit_interior(p)))

# From uniform draws, exact zeros replaced (see `_nonzero_draw`):
@inline rand_impl(ctx::GenContext, ::StdLogistic) = stdconvert(StdLogistic, StdUniform, _nonzero_draw(rand(ctx)))
@inline function batched_rand_impl(ctx::GenContext, ::StdLogistic, sz::SizeLike)
    stdconvert(StdLogistic, StdUniform, _nonzero_draws(_rand_bulk(ctx, sz)))
end

smf(::StdLogistic, x) = logistic(x)
smf(::StdLogistic) = logistic
csmf(::StdLogistic, x) = logistic(-x)
invcsmf(μ::StdLogistic, p) = -invsmf(μ, p)
logsmf(::StdLogistic, x) = -log1pexp(-x)
logcsmf(::StdLogistic, x) = -log1pexp(x)
invlogsmf(::StdLogistic, logp) = (lp = min(logp, zero(logp)); _nan_above_zero(logp, lp - log1mexp(lp)))
invlogcsmf(μ::StdLogistic, logp) = -invlogsmf(μ, logp)

invsmf(μ::StdLogistic, p) = transport_def(μ, StdUniform(), p)

# Interval masses and their inverses in factored forms without
# cancellation, `logistic(x) - logistic(a) = -expm1(a - x) logistic(x)
# logistic(-a)` and `x - a = log1p(p / logistic(a)) - log1p(-p /
# logistic(-a))`, the inverses floored at the quantiles of the smallest
# positive float:
function smf(::StdLogistic, a, x)
    d = ifelse(a == x, zero(a - x), a - x)
    -expm1(d) * logistic(x) * logistic(-a)
end

function invsmf(::StdLogistic, a, p)
    L_a, L_na = logistic(a), logistic(-a)
    valid = (-L_a <= p) & (p <= L_na)
    # From the quantile of `smf(a) + p` for infinite `a`, arguments
    # sanitized for invalid `p` and infinite `a`:
    finite = isfinite(a)
    a_s, p_s = ifelse(finite, a, zero(a)), ifelse(valid & finite, p, zero(p))
    x_a = a_s + (log1p(p_s / logistic(a_s)) - log1p(-p_s / logistic(-a_s)))
    x = ifelse(finite, x_a, logit(clamp(L_a + p, zero(L_a), one(L_a))))
    x_max = -log(_prob_floor(x))
    ifelse(valid, clamp(x, -x_max, x_max), oftype(x, NaN))
end
