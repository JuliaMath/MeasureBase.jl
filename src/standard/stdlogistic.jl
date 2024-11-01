"""
    StdLogistic <: StdMeasure

Represents the standard (centered, scale of one)
[logistic](https://en.wikipedia.org/wiki/Logistic_distribution) probability measure.

See [`StdMeasure`](@ref) for the semantics of standard measures in the
context of MeasureBase.
"""
struct StdLogistic <: StdMeasure end
export StdLogistic

@inline insupport(::StdLogistic, x) = true

@inline logdensity_def(::StdLogistic, x) = (u = -abs(x); u - 2 * log1pexp(u))
@inline basemeasure(::StdLogistic) = LebesgueBase()

@inline transport_def(::StdUniform, μ::StdLogistic, x) = logistic(x)
@inline transport_def(::StdLogistic, μ::StdUniform, p) = logit(p)

@inline function Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdLogistic) where {T}
    logit(rand(rng, T))
end

smf(::StdLogistic, x) = logistic(x)
smf(::StdLogistic) = logistic

invsmf(::StdLogistic, p) = logit(p)
invsmf(::StdLogistic) = logit
