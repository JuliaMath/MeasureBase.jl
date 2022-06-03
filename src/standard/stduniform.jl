struct StandardUniform <: AbstractMeasure end

export StandardUniform

insupport(d::StandardUniform, x) = 0 ≤ x ≤ 1

@inline logdensity_def(::StandardUniform, x) = zero(x)
@inline basemeasure(::StandardUniform) = Lebesgue()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StandardUniform) where {T} = randn(rng, T)

struct StandardUniform <: AbstractMeasure end

export StandardUniform

insupport(d::StandardUniform, x) = zero(x) ≤ x ≤ one(x)

@inline logdensity_def(::StandardUniform, x) = -x^2 / 2
@inline basemeasure(::StandardUniform) = LebesgueMeasure()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StandardUniform) where {T} = rand(rng, T)
