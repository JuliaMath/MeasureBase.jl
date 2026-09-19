
export Dirac

struct Dirac{X} <: AbstractMeasure
    x::X
end

function Pretty.tile(d::Dirac)
    Pretty.literal("Dirac(") * Pretty.tile(d.x) * Pretty.literal(")")
end

Base.:(==)(a::Dirac, b::Dirac) = a.x == b.x
Base.hash(a::Dirac, h::UInt) = hash(a.x, hash(:Dirac, h))
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

# Batched kernels cover Dirac measures with numbers and numeric arrays as
# flat variates, others have no declared variate rank:
const _FlatDirac = Dirac{<:Union{Number,AbstractArray{<:Number}}}

@inline batched_transport_to_std(::Type{S}, ::Dirac, ::Number) where {S<:StdMeasure} = SVector{0,Bool}()
function batched_transport_to_std(::Type{S}, μ::_FlatDirac, X::AbstractArray) where {S<:StdMeasure}
    n = length(_value_flatsize(μ.x))
    similar(X, Bool, (0, ntuple(i -> size(X, n + i), Val(ndims(X) - n))...))
end

function batched_transport_from_std(::Type{S}, μ::_FlatDirac, Z::AbstractArray) where {S<:StdMeasure}
    _const_variates(μ.x, Z)
end
@inline _const_variates(x::Number, ::AbstractVector) = x
function _const_variates(x, Z::AbstractArray)
    X = similar(Z, eltype(x), (size(x)..., Base.tail(size(Z))...))
    X .= x
    return X
end

@inline mspace_ndims(::Type{<:Dirac{<:AbstractArray{<:Number,N}}}) where {N} = N

# Batches of array variates: all elements of a variate must match.
function batched_logdensityof_impl(μ::Dirac{<:AbstractArray{<:Number,N}}, X::AbstractArray) where {N}
    matches = _all_leading_dims(X .== μ.x, static(N))
    ifelse.(matches, zero(_logd_numtype(X)), _neg_inf_logd(X))
end
function batched_logdensity_def(μ::Dirac{<:AbstractArray{<:Number}}, X::AbstractArray)
    _zero_logd_batch(X, static(ndims(μ.x)))
end

@inline _all_leading_dims(A::AbstractArray{Bool,N}, ::StaticInteger{N}) where {N} = all(A)
@inline function _all_leading_dims(A::AbstractArray{Bool}, ::StaticInteger{N}) where {N}
    _drop_leading_dims(all(A; dims = ntuple(identity, Val(N))), static(N))
end

Adapt.adapt_structure(to, μ::Dirac) = Dirac(Adapt.adapt(to, μ.x))
