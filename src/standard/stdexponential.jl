struct StdExponential <: StdMeasure end

export StdExponential

insupport(::StdExponential, x) = x ≥ zero(x)

@inline function logdensityof(::StdExponential, x)
    R = float(typeof(x))
    x ≥ zero(R) ? convert(R, -x) : R(-Inf)
end

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = LebesgueBase()

@inline transport_def(::StdUniform, μ::StdExponential, x) = -expm1(-x)
@inline transport_def(::StdExponential, μ::StdUniform, x) = -log1p(-x)

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdExponential) where {T} = randexp(rng, T)
