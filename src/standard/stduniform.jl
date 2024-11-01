"""
    StdUniform <: StdMeasure

Represents the standard
[uniform](https://en.wikipedia.org/wiki/Continuous_uniform_distribution)
probability measure (from zero to one). It is the
same as the Lebesgue measure restricted to the unit interval.

See [`StdMeasure`](@ref) for the semantics of standard measures in the
context of MeasureBase.
"""
struct StdUniform <: StdMeasure end
export StdUniform

insupport(::StdUniform, x) = zero(x) ≤ x ≤ one(x)

@inline logdensity_def(::StdUniform, x) = zero(x)
@inline basemeasure(::StdUniform) = LebesgueBase()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdUniform) where {T} = rand(rng, T)

massof(::StdUniform, s::Interval) = massof(Lebesgue(𝕀), s::Interval)

smf(::StdUniform, x) = clamp(x, zero(x), one(x))

function invsmf(::StdUniform, p)
    @assert zero(p) ≤ p ≤ one(p)
    p
end
