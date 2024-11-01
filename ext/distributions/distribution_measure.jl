# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).


const DistributionMeasure{F<:VariateForm,S<:ValueSupport,D<:Distribution{F,S}} = AsMeasure{D}

@inline MeasureBase.AbstractMeasure(obj::Distribution) = AsMeasure{typeof(obj)}(obj)
@inline Base.convert(::Type{AbstractMeasure}, obj::Distribution) = AbstractMeasure(obj)

@inline Distributions.Distribution(m::DistributionMeasure) = m.obj
@inline Distributions.Distribution{F}(m::DistributionMeasure{F}) where {F<:VariateForm} = Distribution(m)
@inline Distributions.Distribution{F,S}(m::DistributionMeasure{F,S}) where {F<:VariateForm,S<:ValueSupport} = Distribution(m)

@inline Base.convert(::Type{Distribution}, m::DistributionMeasure) = Distribution(m)
@inline Base.convert(::Type{Distribution{F}}, m::DistributionMeasure{F}) where {F<:VariateForm} = Distribution(m)
@inline Base.convert(::Type{Distribution{F,S}}, m::DistributionMeasure{F,S}) where {F<:VariateForm,S<:ValueSupport} = Distribution(m)


Base.rand(rng::AbstractRNG, ::Type{T}, m::DistributionMeasure) where {T<:Real} = convert_realtype(T, rand(m.obj))

function _flat_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution{<:ArrayLikeVariate{0}}, sz::Dims) where {T<:Real}
    convert_realtype(T, reshape(rand(d, prod(sz)), sz...))
end

function _flat_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution{<:ArrayLikeVariate{1}}, sz::Dims) where {T<:Real}
    convert_realtype(T, reshape(rand(rng, d, prod(sz)), size(d)..., sz...))
end

function _flat_powrand(rng::AbstractRNG, ::Type{T}, d::ReshapedDistribution{N,<:Any,<:Distribution{<:ArrayLikeVariate{1}}}, sz::Dims) where {T<:Real,N}
    convert_realtype(T, reshape(rand(rng, d.dist, prod(sz)), d.dims..., sz...))
end

function _flat_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution, sz::Dims) where {T<:Real}
    flatview(ArrayOfSimilarArrays(convert_realtype(T, rand(rng, d, sz))))
end

function Base.rand(rng::AbstractRNG, ::Type{T}, m::PowerMeasure{<:DistributionMeasure{<:ArrayLikeVariate{0}}, NTuple{N,Base.OneTo{Int}}}) where {T<:Real,N}
    _flat_powrand(rng, T, m.parent.obj, map(length, m.axes))
end

function Base.rand(rng::AbstractRNG, ::Type{T}, m::PowerMeasure{<:DistributionMeasure{<:ArrayLikeVariate{M}}, NTuple{N,Base.OneTo{Int}}}) where {T<:Real,M,N}
    flat_data = _flat_powrand(rng, T, m.parent.obj, map(length, m.axes))
    ArrayOfSimilarArrays{T,M,N}(flat_data)
end


@inline DensityInterface.densityof(m::DistributionMeasure) = densityof(m.obj)
@inline DensityInterface.logdensityof(m::DistributionMeasure) = logdensityof(m.obj)

@inline MeasureBase.logdensity_def(m::DistributionMeasure, x) = DensityInterface.logdensityof(m.obj, x)
@inline MeasureBase.unsafe_logdensityof(m::DistributionMeasure, x) = DensityInterface.logdensityof(m.obj, x)
@inline MeasureBase.insupport(m::DistributionMeasure, x) = Distributions.insupport(m.obj, x)

@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate{0},<:Continuous}) = Lebesgue()
@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate,<:Continuous}) = Lebesgue()^size(m.obj)
@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate{0},<:Discrete}) = Counting()
@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate,<:Discrete}) = Counting()^size(m.obj)

@inline MeasureBase.basemeasure(m::DistributionMeasure) = rootmeasure(m)

@inline MeasureBase.mspace_elsize(m::DistributionMeasure{<:ArrayLikeVariate}) = size(m.obj)

@inline MeasureBase.getdof(m::DistributionMeasure{<:ArrayLikeVariate{0}}) = 1

@inline MeasureBase.paramnames(m::DistributionMeasure) = propertynames(m.obj)
@inline MeasureBase.params(m::DistributionMeasure) = NamedTuple{propertynames(m.obj)}(Distributions.params(m.obj))

# @inline MeasureBase.testvalue(m::DistributionMeasure) = testvalue(basemeasure(d))


@inline MeasureBase.basemeasure(d::Distributions.Poisson) = Counting(MeasureBase.BoundedInts(static(0), static(Inf)))
@inline MeasureBase.basemeasure(d::Distributions.Product{<:Any,<:Distributions.Poisson}) = Counting(MeasureBase.BoundedInts(static(0), static(Inf)))^size(d)
