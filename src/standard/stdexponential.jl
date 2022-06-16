struct StdExponential <: StdMeasure end

export StdExponential

insupport(d::StdExponential, x) = x ≥ zero(x)

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = Lebesgue()

@inline getdof(::StdExponential) = static(1)

function Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdExponential) where {T}
    randexp(rng, T)
end
