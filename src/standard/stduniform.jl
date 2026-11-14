"""
    StdUniform <: StdMeasure

The standard uniform measure on the unit interval, the uniform
distribution on `[0, 1]` as a measure.
"""
struct StdUniform <: StdMeasure end

export StdUniform

insupportof(::StdUniform, x) = (zero(x) ≤ x) & (x ≤ one(x))

@inline function logdensityof_impl(d::StdUniform, x)
    R = float(typeof(x))
    _checksupport(insupportof(d, x), zero(R))
end

@inline logdensity_def(::StdUniform, x) = zero(x)
@inline basemeasure(::StdUniform) = LebesgueBase()

@inline rand_impl(ctx::GenContext, ::StdUniform) = rand(ctx)
@inline batched_rand_impl(ctx::GenContext, ::StdUniform, sz::SizeLike) = _rand_bulk(ctx, sz)

massof(::StdUniform, s::Interval) = massof(Lebesgue(0.0 .. 1.0), s)

smf(::StdUniform, x) = clamp(x, zero(x), one(x))
csmf(::StdUniform, x) = clamp(one(x) - x, zero(x), one(x))
invcsmf(d::StdUniform, p) = _nan_outside(d, p, one(p) - p)

invsmf(d::StdUniform, p) = _nan_outside(d, p, p)

# A draw of exact zero stands for the lowest cell of the RNG grid, not for
# its lower bound, which `_unit_interior` would clamp to `floatmin` (an
# atom in the transported variates). Uniform draws for transports replace
# exact zeros by half the smallest nonzero draw of the RNG grid of `rand`
# (`eps(T) / 4`), other draws stay as they are:
@inline _nonzero_draw(u::Number) = ifelse(iszero(u), eps(typeof(one(u))) / 4, u)
@inline _nonzero_draws(u::Number) = _nonzero_draw(u)
@inline _nonzero_draws(U::AbstractArray) = _nonzero_draw.(U)

@inline _rand_std(ctx::GenContext, ::Type{StdUniform}, dims::SizeLike) = _nonzero_draws(batched_rand_impl(ctx, StdUniform(), dims))
