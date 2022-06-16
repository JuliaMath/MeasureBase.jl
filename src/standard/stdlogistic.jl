struct StdLogistic <: StdMeasure end

export StdLogistic

@inline insupport(d::StdLogistic, x) = true

@inline logdensity_def(::StdLogistic, x) = (u = -abs(x); u - 2*log1pexp(u))
@inline basemeasure(::StdLogistic) = Lebesgue()

@inline getdof(::StdLogistic) = static(1)

@inline vartransform_def(::StdUniform, μ::StdLogistic, x) = logistic(checked_var(μ, x))
@inline vartransform_def(::StdLogistic, μ::StdUniform, x) = logit(checked_var(μ, x))

@inline Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdLogistic) where {T} = logit(rand(rng, T))
