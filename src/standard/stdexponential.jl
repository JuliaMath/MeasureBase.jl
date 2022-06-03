struct StandardExponential <: AbstractMeasure end

export StandardExponential

insupport(d::StandardExponential, x) = x ≥ zero(x)

@inline logdensity_def(::StandardExponential, x) = -x
@inline basemeasure(::StandardExponential) = Lebesgue()

function Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StandardExponential) where {T}
    randexp(rng, T)
end
