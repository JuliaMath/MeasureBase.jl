
"""
    proxy(μ::AbstractMeasure)::AbstractMeasure

Return a measure equivalent to `μ` to which methods like `logdensity_def`
and `basemeasure` of `μ` can be delegated, via [`MeasureBase.@useproxy`](@ref).

For example, a `Normal{(:μ,:σ)}` is equivalent to an affine transformation
of a `Normal{()}`. Using a proxy instead of constructing the transformed
measure directly keeps dispatch on the original measure type.
"""
function proxy end

"""
    MeasureBase.@useproxy M

Delegate the measure interface for measures of type `M` to their
[`proxy`](@ref).
"""
macro useproxy(M)
    M = esc(M)
    quote
        @inline $MeasureBase.logdensity_def(μ::$M, x) = logdensity_def(proxy(μ), x)
        @inline $MeasureBase.unsafe_logdensityof(μ::$M, x) = unsafe_logdensityof(proxy(μ), x)

        @inline $MeasureBase.basemeasure(μ::$M) = basemeasure(proxy(μ))
        @inline $MeasureBase.basemeasure_depth(μ::$M) = basemeasure_depth(proxy(μ))
        @inline $MeasureBase.rootmeasure(μ::$M) = rootmeasure(proxy(μ))

        @inline $MeasureBase.insupportof(μ::$M) = insupportof(proxy(μ))
        @inline $MeasureBase.insupportof(μ::$M, x) = insupportof(proxy(μ), x)

        @inline $MeasureBase.getdof(μ::$M) = getdof(proxy(μ))
        @inline $MeasureBase.fast_dof(μ::$M) = fast_dof(proxy(μ))

        @inline $MeasureBase.localmeasure(μ::$M, x) = localmeasure(proxy(μ), x)
        @inline $MeasureBase.transportmeasure(μ::$M, x) = transportmeasure(proxy(μ), x)

        @inline $MeasureBase.massof(μ::$M) = massof(proxy(μ))
        @inline $MeasureBase.massof(μ::$M, s) = massof(proxy(μ), s)
        # Disambiguation with massof(μ, ::AbstractInterval):
        @inline $MeasureBase.massof(μ::$M, s::$(IntervalSets.AbstractInterval)) =
            massof(proxy(μ), s)

        @inline $MeasureBase.smf(μ::$M, x) = smf(proxy(μ), x)
        @inline $MeasureBase.invsmf(μ::$M, x) = invsmf(proxy(μ), x)
        (μ::$M)(s) = proxy(μ)(s)
    end
end
