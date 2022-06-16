struct StdLogistic <: StdMeasure end

export StdLogistic

@inline _isreal(x) = false
@inline _isreal(x::Real) = true

@inline insupport(d::StdLogistic, x) = _isreal(x)
@inline insupport(d::StdLogistic) = _isreal

@inline logdensity_def(::StdLogistic, x) = (u = -abs(x); u - 2*log1pexp(u))
@inline basemeasure(::StdLogistic) = Lebesgue()

@inline getdof(::StdLogistic) = static(1)

@inline vartransform_def(::StdUniform, ::StdLogistic, x::Real) = logistic(x)
@inline vartransform_def(::StdLogistic, ::StdUniform, x::Real) = logit(x)

@inline Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdLogistic) where {T} = logit(rand(rng, T))

@inline StdMeasure(::typeof(randn)) = StdLogistic()
