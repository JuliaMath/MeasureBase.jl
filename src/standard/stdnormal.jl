struct StandardNormal <: AbstractMeasure end

export StandardNormal

insupport(d::StandardNormal, x) = true
insupport(d::StandardNormal) = Returns(true)

@inline logdensity_def(::StandardNormal, x) = -x^2 / 2
@inline basemeasure(::StandardNormal) = WeightedMeasure(static(-0.5 * log2π), Lebesgue(ℝ))

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StandardNormal) where {T} = randn(rng, T)
