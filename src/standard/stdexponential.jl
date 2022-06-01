struct StdExponential <: AbstractMeasure end

export StdExponential

insupport(d::StdExponential, x) = x ≥ zero(x)

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = Lebesgue()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdExponential) where {T} = randexp(rng, T)

