struct StdExponential <: StdMeasure end

export StdExponential

insupport(d::StdExponential, x) = x ≥ zero(x)

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = LebesgueBase()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdExponential) where {T} = randexp(rng, T)

transport_origin(::StdExponential) = StdUniform()

from_origin(::StdExponential, x) = -log1p(-x)
to_origin(::StdExponential, x) = -expm1(-x)