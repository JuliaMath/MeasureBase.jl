struct StdExponential <: AbstractMeasure end

export StdExponential

insupport(d::StdExponential, x) = x ≥ zero(x)

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = Lebesgue()

function Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdExponential) where {T}
    randexp(rng, T)
end
