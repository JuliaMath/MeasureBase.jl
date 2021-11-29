
"""
function proxy end

It's often useful to delegate methods like `logdensity` and `basemeasure` to
those of a different measure. For example, a `Normal{(:μ,:σ)}` is equivalent to
an affine transformation of a `Normal{()}`.

We _could_ just have calls like `Normal(μ=2,σ=4)` directly construct a
transformed measure, but this would make dispatch awkward.
"""
function proxy end

proxy(μ) = μ

proxy(f, μ) = proxy(μ)

logdensity_def(μ, x) = logdensity_def(proxy(μ), x)


basemeasure(μ, x) = basemeasure(proxy(μ), x)

basemeasure(μ) = basemeasure(proxy(μ))
