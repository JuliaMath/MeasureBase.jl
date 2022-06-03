struct StdUniform <: AbstractMeasure end

export StdUniform

insupport(d::StdUniform, x) = 0 ≤ x ≤ 1

@inline logdensity_def(::StdUniform, x) = zero(x)
@inline basemeasure(::StdUniform) = Lebesgue()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdUniform) where {T} = randn(rng, T)
