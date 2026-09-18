
export Dirac

struct Dirac{X} <: AbstractMeasure
    x::X
end

function Pretty.tile(d::Dirac)
    Pretty.literal("Dirac(") * Pretty.tile(d.x) * Pretty.literal(")")
end

Base.:(==)(a::Dirac, b::Dirac) = a.x == b.x
Base.isapprox(a::Dirac, b::Dirac; kwargs...) = isapprox(a.x, b.x; kwargs...)

gentype(μ::Dirac{X}) where {X} = X

function (μ::Dirac{X})(s) where {X}
    μ.x ∈ s && return 1
    return 0
end

basemeasure(d::Dirac) = CountingBase()

massof(::Dirac) = static(1.0)

function logdensityof_impl(μ::Dirac, x::Number)
    R = float(typeof(x))
    _checksupport(insupport(μ, x), zero(R))
end

logdensityof_impl(μ::Dirac, x) = _checksupport(insupport(μ, x), zero(_logd_numtype(x)))

logdensity_def(::Dirac, x::Number) = zero(float(typeof(x)))
logdensity_def(::Dirac, x) = zero(_logd_numtype(x))

@inline rand_impl(::GenContext, μ::Dirac) = μ.x
@inline batched_rand_impl(ctx::GenContext, μ::Dirac, sz::Dims) = _const_batch(ctx, μ.x, sz)

export dirac

dirac(d::AbstractMeasure) = Dirac(rand(d))

insupport(d::Dirac, x) = x == d.x

@inline getdof(::Dirac) = static(0)

@inline mspace_elsize(μ::Dirac) = _value_elsize(μ.x)
@inline mspace_flatsize(μ::Dirac) = _value_flatsize(μ.x)
@inline mspace_flatsize(::Type{<:Dirac{<:Number}}) = ()

@propagate_inbounds function checked_arg(μ::Dirac, x)
    @boundscheck insupport(μ, x) || throw(ArgumentError("Invalid variate for measure"))
    x
end

# Dirac measures have no degrees of freedom:
@inline transport_to_std(::Type{S}, ::Dirac, x) where {S<:StdMeasure} = SVector{0,Bool}()
@inline transport_from_std(::Type{S}, μ::Dirac, z::AbstractVector) where {S<:StdMeasure} = μ.x
@inline transport_from_std_with_rest(::Type{S}, μ::Dirac, z::AbstractVector) where {S<:StdMeasure} = μ.x, z

function batched_transport_to_std(::Type{S}, μ::Dirac, X::AbstractArray) where {S<:StdMeasure}
    n = length(_value_flatsize(μ.x))
    similar(X, Bool, (0, ntuple(i -> size(X, n + i), Val(ndims(X) - n))...))
end

function batched_transport_from_std(::Type{S}, μ::Dirac, Z::AbstractArray) where {S<:StdMeasure}
    X = similar(Z, eltype(μ.x), (size(μ.x)..., Base.tail(size(Z))...))
    X .= μ.x
    return X
end

@inline mspace_ndims(::Type{<:Dirac{<:AbstractArray{<:Any,N}}}) where {N} = N

# Batches of array variates: all elements of a variate must match.
function batched_logdensityof_impl(μ::Dirac{<:AbstractArray{<:Any,N}}, X::AbstractArray) where {N}
    matches = _all_leading_dims(X .== μ.x, static(N))
    ifelse.(matches, zero(_logd_numtype(X)), _neg_inf_logd(X))
end
batched_logdensity_def(μ::Dirac{<:AbstractArray}, X::AbstractArray) = _zero_logd_batch(X, static(ndims(μ.x)))

@inline _all_leading_dims(A::AbstractArray{Bool,N}, ::StaticInteger{N}) where {N} = all(A)
@inline function _all_leading_dims(A::AbstractArray{Bool}, ::StaticInteger{N}) where {N}
    dropdims(all(A; dims = ntuple(identity, Val(N))); dims = ntuple(identity, Val(N)))
end

Adapt.adapt_structure(to, μ::Dirac) = Dirac(Adapt.adapt(to, μ.x))
