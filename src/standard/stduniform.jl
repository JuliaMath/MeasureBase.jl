struct StdUniform <: StdMeasure end

export StdUniform

insupport(d::StdUniform, x) = zero(x) ≤ x ≤ one(x)

@inline logdensity_def(::StdUniform, x) = zero(x)
@inline basemeasure(::StdUniform) = Lebesgue()

@inline getdof(::StdUniform) = static(1)

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdUniform) where {T} = randn(rng, T)
