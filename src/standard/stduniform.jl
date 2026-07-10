struct StdUniform <: StdMeasure end

export StdUniform

insupport(::StdUniform, x) = (zero(x) ≤ x) & (x ≤ one(x))

@inline function logdensityof_impl(d::StdUniform, x)
    R = float(typeof(x))
    _checksupport(insupport(d, x), zero(R))
end

@inline logdensity_def(::StdUniform, x) = zero(x)
@inline basemeasure(::StdUniform) = LebesgueBase()

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdUniform) where {T} = rand(rng, T)

massof(::StdUniform, s::Interval) = massof(Lebesgue(0.0 .. 1.0), s)

smf(::StdUniform, x) = clamp(x, zero(x), one(x))

function invsmf(::StdUniform, p)
    @assert zero(p) ≤ p ≤ one(p)
    p
end
