struct StdNormal <: AbstractMeasure end

export StdNormal

insupport(d::StdNormal, x) = true
insupport(d::StdNormal) = Returns(true)

@inline logdensity_def(::StdNormal, x) = -x^2 / 2
@inline basemeasure(::StdNormal) = WeightedMeasure(static(-0.5 * log2π), Lebesgue(ℝ))

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdNormal) where {T} = randn(rng, T)
