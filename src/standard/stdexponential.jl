struct StdExponential <: StdMeasure end

export StdExponential

insupport(d::StdExponential, x) = x ≥ zero(x)

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = Lebesgue()

@inline getdof(::StdExponential) = static(1)

@inline vartransform_def(::StdUniform, μ::StdExponential, x) = - expm1(- checked_var(μ, x))
@inline vartransform_def(::StdExponential, μ::StdUniform, x) = - log1p(- checked_var(μ, x))

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdExponential) where {T} = randexp(rng, T)

