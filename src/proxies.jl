
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
        #!!!!!!!!!!! TODO add new API methods like localmeasure, transportmeasure, etc. !!!!!!!!!!!!!

        @inline $MeasureBase.logdensity_def(μ::$M, x) = logdensity_def(proxy(μ), x)

        @inline $MeasureBase.basemeasure(μ::$M) = basemeasure(proxy(μ))

        @inline $MeasureBase.basemeasure_depth(μ::$M) = basemeasure_depth(proxy(μ))

        @inline $MeasureBase.transport_origin(μ::$M) = transport_origin(proxy(μ))
        @inline $MeasureBase.to_origin(μ::$M, y) = to_origin(proxy(μ), y)
        @inline $MeasureBase.from_origin(μ::$M, x) = from_origin(proxy(μ), x)

        @inline $MeasureBase.massof(μ::$M) = massof(proxy(μ))
        @inline $MeasureBase.massof(μ::$M, s) = massof(proxy(μ), s)

        @inline $MeasureBase.smf(μ::$M, x) = smf(proxy(μ), x)
        @inline $MeasureBase.invsmf(μ::$M, x) = invsmf(proxy(μ), x)
        (μ::$M)(s) = proxy(μ)(s)
    end
end
