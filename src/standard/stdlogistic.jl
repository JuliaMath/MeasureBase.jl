struct StdLogistic <: StdMeasure end

export StdLogistic

@inline insupport(d::StdLogistic, x) = true

@inline logdensity_def(::StdLogistic, x) = (u = -abs(x); u - 2 * log1pexp(u))
@inline basemeasure(::StdLogistic) = LebesgueBase()

@inline transport_def(::StdUniform, μ::StdLogistic, x) = smf(StdLogistic(), x)
@inline transport_def(::StdLogistic, μ::StdUniform, p) = invsmf(StdLogistic(), p)

@inline function Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdLogistic) where {T}
    logit(rand(rng, T))
end

smf(::StdLogistic, x) = logistic(x)

invsmf(::StdLogistic, p) = logit(p)
