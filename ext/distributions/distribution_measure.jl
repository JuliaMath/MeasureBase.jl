# This file is a part of DistributionMeasures.jl, licensed under the MIT License (MIT).


"""
    struct DistributionMeasure <: AbstractMeasure

Wraps a `Distributions.Distribution` as a `MeasureBase.AbstractMeasure`.

Avoid calling `DistributionMeasure(d::Distribution)` directly. Instead, use
`AbstractMeasure(d::Distribution)` to allow for specialized `Distribution`
to `AbstractMeasure` conversions.

Use `convert(Distribution, m::DistributionMeasure)` or
`Distribution(m::DistributionMeasure)` to convert back to a `Distribution`.
"""
struct DistributionMeasure{F<:VariateForm,S<:ValueSupport,D<:Distribution{F,S}} <: AbstractMeasure
    d::D
end


@inline MeasureBase.AbstractMeasure(d::Distribution) = DistributionMeasure(d)

@inline Base.convert(::Type{AbstractMeasure}, d::Distribution) = DistributionMeasure(d)

@inline Distributions.Distribution(m::DistributionMeasure) = m.d
@inline Distributions.Distribution{F}(m::DistributionMeasure{F}) where {F<:VariateForm} = Distribution(m)
@inline Distributions.Distribution{F,S}(m::DistributionMeasure{F,S}) where {F<:VariateForm,S<:ValueSupport} = Distribution(m)

@inline Base.convert(::Type{Distribution}, m::DistributionMeasure) = Distribution(m)
@inline Base.convert(::Type{Distribution{F}}, m::DistributionMeasure{F}) where {F<:VariateForm} = Distribution(m)
@inline Base.convert(::Type{Distribution{F,S}}, m::DistributionMeasure{F,S}) where {F<:VariateForm,S<:ValueSupport} = Distribution(m)


@inline DensityInterface.densityof(μ::DistributionMeasure) = DensityInterface.densityof(μ.d)
@inline DensityInterface.densityof(μ::DistributionMeasure, x) = DensityInterface.densityof(μ.d, x)
@inline DensityInterface.logdensityof(μ::DistributionMeasure) = DensityInterface.logdensityof(μ.d)
@inline DensityInterface.logdensityof(μ::DistributionMeasure, x) = DensityInterface.logdensityof(μ.d, x)

@inline MeasureBase.logdensity_def(μ::DistributionMeasure, x) = MeasureBase.logdensity_def(μ.d, x)
@inline MeasureBase.unsafe_logdensityof(μ::DistributionMeasure, x) = MeasureBase.unsafe_logdensityof(μ.d, x)
@inline MeasureBase.insupport(μ::DistributionMeasure, x) = MeasureBase.insupport(μ.d, x)
@inline MeasureBase.basemeasure(μ::DistributionMeasure) = MeasureBase.basemeasure(μ.d)
@inline MeasureBase.paramnames(μ::DistributionMeasure) = MeasureBase.paramnames(μ.d)
@inline MeasureBase.params(μ::DistributionMeasure) = MeasureBase.params(μ.d)
@inline MeasureBase.transport_origin(ν::DistributionMeasure) = ν.d
@inline MeasureBase.to_origin(::DistributionMeasure, y) = y
@inline MeasureBase.from_origin(::DistributionMeasure, x) = x


Base.rand(rng::AbstractRNG, ::Type{T}, m::DistributionMeasure) where {T<:Real} = convert_realtype(T, rand(m.d))

function _flat_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution{<:ArrayLikeVariate{0}}, sz::Dims) where {T<:Real}
    convert_realtype(T, reshape(rand(d, prod(sz)), sz...))
end

function _flat_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution{<:ArrayLikeVariate{1}}, sz::Dims) where {T<:Real}
    convert_realtype(T, reshape(rand(d, prod(sz)), size(d)..., sz...))
end

function _flat_powrand(rng::AbstractRNG, ::Type{T}, d::ReshapedDistribution{N,<:Any,<:Distribution{<:ArrayLikeVariate{1}}}, sz::Dims) where {T<:Real,N}
    convert_realtype(T, reshape(rand(d.dist, prod(sz)), d.dims..., sz...))
end

function _flat_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution, sz::Dims) where {T<:Real}
    flatview(ArrayOfSimilarArrays(convert_realtype(T, rand(d, sz))))
end

function Base.rand(rng::AbstractRNG, ::Type{T}, m::PowerMeasure{<:DistributionMeasure{<:ArrayLikeVariate{0}}, NTuple{N,Base.OneTo{Int}}}) where {T<:Real,N}
    _flat_powrand(rng, T, m.parent.d, map(length, m.axes))
end

function Base.rand(rng::AbstractRNG, ::Type{T}, m::PowerMeasure{<:DistributionMeasure{<:ArrayLikeVariate{M}}, NTuple{N,Base.OneTo{Int}}}) where {T<:Real,M,N}
    flat_data = _flat_powrand(rng, T, m.parent.d, map(length, m.axes))
    ArrayOfSimilarArrays{T,M,N}(flat_data)
end
