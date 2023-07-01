"""
    StdExponential <: StdMeasure

Represents the standard (rate of one)
[exponential](https://en.wikipedia.org/wiki/Exponential_distribution) probability measure.

See [`StdMeasure`](@ref) for the semantics of standard measures in the
context of MeasureBase.
"""
struct StdExponential <: StdMeasure end

export StdExponential

insupport(::StdExponential, x) = x ≥ zero(x)

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = LebesgueBase()

@inline transport_def(::StdUniform, μ::StdExponential, x) = -expm1(-x)
@inline transport_def(::StdExponential, μ::StdUniform, x) = -log1p(-x)

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdExponential) where {T} = randexp(rng, T)
