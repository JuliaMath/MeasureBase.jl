
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
        @inline $__module__.logdensity_def(μ::$M, x) = logdensity_def(proxy(μ), x)

        @inline function $__module__.basemeasure(μ::$M)
            p = proxy(μ)
            b = basemeasure(p)
            return b
        end

        @inline $__module__.basemeasure_depth(μ::$M) = basemeasure_depth(proxy(μ))
    end
end
