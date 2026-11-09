# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Mixtures become superpositions, except continuous univariate mixtures,
# which keep the cdf/quantile-based transport of wrapped distributions:
function MeasureBase.AbstractMeasure(d::Distributions.AbstractMixtureModel)
    superpose(map((w, c) -> w * asmeasure(c), Distributions.probs(d), Distributions.components(d)))
end
MeasureBase.AbstractMeasure(d::Distributions.AbstractMixtureModel{Univariate,Continuous}) = AsMeasure{typeof(d)}(d)


# The tails of continuous univariate mixtures are sums over the tails of
# the components, in log space (Distributions' complements cancel), their
# quantiles the roots of these, between the quantiles of the components:
@inline _trafo_logcdf(d::_ContinuousMixture, x::Number) = _mixture_logtail(_trafo_logcdf, d, x)
@inline _trafo_logccdf(d::_ContinuousMixture, x::Number) = _mixture_logtail(_trafo_logccdf, d, x)
function _mixture_logtail(f::F, d::_ContinuousMixture, x) where {F}
    logsumexp(map((w, c) -> log(w) + f(c, x), Distributions.probs(d), Distributions.components(d)))
end

_dist_quantile(d::_ContinuousMixture, p::Real) = _mixture_quantile(_trafo_logcdf, _dist_quantile, d, p)
_dist_cquantile(d::_ContinuousMixture, q::Real) = _mixture_quantile(_trafo_logccdf, _dist_cquantile, d, q)

# The bisection isn't differentiated, a Newton step on the log tail from
# its result gives the derivatives of the root (by the implicit function
# theorem), with respect to the probability and the mixture parameters:
function _mixture_quantile(logtail::F, quantile::G, d::_ContinuousMixture, p::Real) where {F,G}
    x0 = @_adignore _mixture_bisection(logtail, quantile, d, p)
    slope = @_adignore (logtail === _trafo_logcdf ? 1 : -1) * exp(Distributions.logpdf(d, x0) - logtail(d, x0))
    isfinite(slope) && !iszero(slope) ? x0 - (logtail(d, x0) - log(p)) / slope : x0
end

function _mixture_bisection(logtail::F, quantile::G, d::_ContinuousMixture, p::Real) where {F,G}
    lo, hi = extrema(c -> quantile(c, p), Distributions.components(d))
    lp = log(p)
    # The tail is monotonic, bisection keeps the bracket:
    for _ in 1:200
        mid = lo / 2 + hi / 2
        (mid == lo || mid == hi) && break
        if (logtail(d, mid) < lp) == (logtail === _trafo_logcdf)
            lo = mid
        else
            hi = mid
        end
    end
    lo / 2 + hi / 2
end


const _MixtureMeasure = SuperpositionMeasure{
    <:Union{Tuple{Vararg{WeightedMeasure}},AbstractVector{<:WeightedMeasure}},
}

_mixture_component(m::AsMeasure{<:Distribution}) = m.obj

function Distributions.Distribution(m::_MixtureMeasure)
    components = map(c -> _mixture_component(c.base), collect(values(m.components)))
    prior = map(c -> exp(c.logweight), collect(values(m.components)))
    Distributions.MixtureModel(components, prior)
end

Base.convert(::Type{Distribution}, m::_MixtureMeasure) = Distributions.Distribution(m)
