"""
    StdExponential <: StdMeasure

The standard exponential measure, the exponential distribution with unit
scale as a measure.
"""
struct StdExponential <: StdMeasure end

export StdExponential

insupport(::StdExponential, x) = x ≥ zero(x)

@inline function logdensityof_impl(d::StdExponential, x)
    R = float(typeof(x))
    _checksupport(insupport(d, x), convert(R, -x))
end

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = LebesgueBase()

@inline transport_def(::StdUniform, μ::StdExponential, x) = _nan_outside(μ, x, -expm1(-x))
@inline transport_def(::StdExponential, μ::StdUniform, x) = _nan_outside(μ, x, -log1p(-_unit_interior(x)))

@inline rand_impl(ctx::GenContext, ::StdExponential) = randexp(get_rng(ctx), get_precision(ctx))
@inline batched_rand_impl(ctx::GenContext, ::StdExponential, sz::SizeLike) = _randexp_bulk(ctx, sz)
