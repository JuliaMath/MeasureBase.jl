
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

Base.rand(::Random.AbstractRNG, T::Type, μ::Dirac) = μ.x

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
