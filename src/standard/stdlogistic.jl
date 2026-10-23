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

# Uniform draws may be zero, the uniform transport keeps the draws finite:
@inline rand_impl(ctx::GenContext, ::StdLogistic) = stdconvert(StdLogistic, StdUniform, rand(ctx))
@inline function batched_rand_impl(ctx::GenContext, ::StdLogistic, sz::SizeLike)
    stdconvert(StdLogistic, StdUniform, _rand_bulk(ctx, sz))
end

smf(::StdLogistic, x) = logistic(x)
smf(::StdLogistic) = logistic

invsmf(μ::StdLogistic, p) = transport_def(μ, StdUniform(), p)
