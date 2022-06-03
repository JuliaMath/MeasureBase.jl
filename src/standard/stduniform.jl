struct StdUniform <: AbstractMeasure end

export StdUniform

insupport(d::StdUniform, x) = zero(x) ≤ x ≤ one(x)

@inline logdensity_def(::StdUniform, x) = zero(x)
@inline basemeasure(::StdUniform) = Lebesgue()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdUniform) where {T} = randn(rng, T)
