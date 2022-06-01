struct StdUniform <: AbstractMeasure end

export StdUniform

insupport(d::StdUniform, x) = 0 ≤ x ≤ 1

@inline logdensity_def(::StdUniform, x) = zero(x)
@inline basemeasure(::StdUniform) = Lebesgue()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdUniform) where {T} = randn(rng, T)

struct StdUniform <: AbstractMeasure end

export StdUniform

insupport(d::StdUniform, x) = zero(x) ≤ x ≤ one(x)

@inline logdensity_def(::StdUniform, x) = -x^2 / 2
@inline basemeasure(::StdUniform) = LebesgueMeasure()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdUniform) where {T} = rand(rng, T)

