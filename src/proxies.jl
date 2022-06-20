
"""
function proxy end

It's often useful to delegate methods like `logdensity` and `basemeasure` to
those of a different measure. For example, a `Normal{(:μ,:σ)}` is equivalent to
an affine transformation of a `Normal{()}`.

We _could_ just have calls like `Normal(μ=2,σ=4)` directly construct a
transformed measure, but this would make dispatch awkward.
"""
function proxy end

macro useproxy(M)
    M = esc(M)
    quote
        @inline $MeasureBase.logdensity_def(μ::$M, x) = logdensity_def(proxy(μ), x)

        @inline $MeasureBase.basemeasure(μ::$M) = basemeasure(proxy(μ))

        @inline $MeasureBase.basemeasure_depth(μ::$M) = basemeasure_depth(proxy(μ))
    end
end
