"""
    StdLogistic <: StdMeasure

The standard logistic measure, the logistic distribution with zero
location and unit scale as a measure.
"""
struct StdLogistic <: StdMeasure end

export StdLogistic

@inline insupport(d::StdLogistic, x) = true

@inline logdensityof_impl(::StdLogistic, x) = (u = -abs(x); u - 2 * log1pexp(u))

@inline logdensity_def(::StdLogistic, x) = logdensityof(StdLogistic(), x)
@inline basemeasure(::StdLogistic) = LebesgueBase()

@inline transport_def(::StdUniform, μ::StdLogistic, x) = logistic(x)
@inline transport_def(::StdLogistic, μ::StdUniform, p) = logit(p)

@inline rand_impl(ctx::GenContext, ::StdLogistic) = logit(rand(get_rng(ctx), get_precision(ctx)))
@inline batched_rand_impl(ctx::GenContext, ::StdLogistic, sz::Dims) = logit.(_rand_bulk(ctx, sz))

smf(::StdLogistic, x) = logistic(x)
smf(::StdLogistic) = logistic

invsmf(::StdLogistic, p) = logit(p)
invsmf(::StdLogistic) = logit
