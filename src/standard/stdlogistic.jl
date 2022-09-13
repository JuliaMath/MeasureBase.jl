struct StdLogistic <: StdMeasure end

export StdLogistic

@inline insupport(d::StdLogistic, x) = true

@inline logdensity_def(::StdLogistic, x) = (u = -abs(x); u - 2 * log1pexp(u))
@inline basemeasure(::StdLogistic) = LebesgueMeasure()

@inline transport_def(::StdUniform, μ::StdLogistic, x) = logistic(x)
@inline transport_def(::StdLogistic, μ::StdUniform, x) = logit(x)

@inline function Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdLogistic) where {T}
    logit(rand(rng, T))
end
