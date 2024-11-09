export Half

struct Half{M} <: AbstractMeasure
    parent::M
end

function Base.show(io::IO, μ::Half)
    print(io, "Half")
    show(io, μ.parent)
end

unhalf(μ::Half) = μ.parent

@inline function basemeasure(μ::Half)
    weightedmeasure(logtwo, basemeasure(unhalf(μ)))
end

function Base.rand(rng::AbstractRNG, ::Type{T}, μ::Half) where {T}
    return abs(rand(rng, T, unhalf(μ)))
end

function strict_logdensityof(μ::Half, x)
    ld = strict_logdensityof(unhalf(μ), x) - loghalf
    return x ≥ 0 ? ld : oftype(ld, -Inf)
end

logdensity_def(μ::Half, x) = logdensity_def(unhalf(μ), x)

@inline function insupport(d::Half, x)
    x ≥ 0 || return false
    insupport(unhalf(d), x)
end

testvalue(::Type{T}, ::Half) where {T} = one(T)

massof(μ::Half) = massof(unhalf(μ))

function smf(μ::Half, x)
    2 * smf(μ.parent, max(x, zero(x))) - 1
end

function invsmf(μ::Half, p)
    @assert zero(p) ≤ p ≤ one(p)
    invsmf(μ.parent, (p + 1) / 2)
end

transport_def(μ::Half, ::StdUniform, p) = invsmf(μ, p)
transport_def(::StdUniform, μ::Half, x) = smf(μ, x)
