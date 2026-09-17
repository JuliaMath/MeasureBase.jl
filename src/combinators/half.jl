export Half

struct Half{M} <: AbstractMeasure
    parent::M
end

@inline mspace_elsize(μ::Half) = mspace_elsize(μ.parent)
@inline mspace_flatsize(μ::Half) = mspace_flatsize(μ.parent)
@inline mspace_flatsize(::Type{<:Half{M}}) where {M} = mspace_flatsize(M)
@inline preferred_stdmeasure(::Type{<:Half}) = StdUniform

function Base.show(io::IO, μ::Half)
    print(io, "Half")
    show(io, μ.parent)
end

unhalf(μ::Half) = μ.parent

@inline function basemeasure(μ::Half)
    weightedmeasure(logtwo, basemeasure(unhalf(μ)))
end

@inline rand_impl(ctx::GenContext, μ::Half) = abs(rand_impl(ctx, unhalf(μ)))
@inline batched_rand_impl(ctx::GenContext, μ::Half, sz::Dims) = abs.(batched_rand_impl(ctx, unhalf(μ), sz))

function logdensityof_impl(μ::Half, x)
    ld = logdensityof(unhalf(μ), x) - loghalf
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

@inline transport_to_std(::Type{StdUniform}, μ::Half, x) = smf(μ, x)
@inline transport_from_std(::Type{StdUniform}, μ::Half, p) = invsmf(μ, p)
