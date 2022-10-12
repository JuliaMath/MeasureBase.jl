struct StdUniform <: StdMeasure end

export StdUniform

insupport(d::StdUniform, x) = zero(x) ≤ x ≤ one(x)

@inline logdensity_def(::StdUniform, x) = zero(x)
@inline basemeasure(::StdUniform) = LebesgueBase()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdUniform) where {T} = rand(rng, T)

massof(::StdUniform, s::Interval) = massof(Lebesgue(𝕀), s::Interval)

smf(::StdUniform, x) = clamp(x, zero(x), one(x))

function smfinv(::StdUniform, p)
    @assert zero(p) ≤ p ≤ one(p)
    p
end
