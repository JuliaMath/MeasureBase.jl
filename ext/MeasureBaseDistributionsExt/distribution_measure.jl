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


# Distributions' samplers run on the CPU, variates on other compute units
# are generated from standard variates via the transports:
MeasureBase.rand_impl(ctx::GenContext, m::DistributionMeasure) = _dist_rand(ctx, m, get_compute_unit(ctx))
MeasureBase.batched_rand_impl(ctx::GenContext, m::DistributionMeasure, sz::Dims) = _dist_batched_rand(ctx, m, sz, get_compute_unit(ctx))

_dist_rand(ctx::GenContext, m::DistributionMeasure, ::CPUnit) =
    convert_realtype(get_precision(ctx), rand(get_rng(ctx), m.obj))
_dist_rand(ctx::GenContext, m::DistributionMeasure, ::AbstractComputeUnit) =
    MeasureBase._rand_default(ctx, m, (), MeasureBase._NoRandImpl())
_dist_batched_rand(ctx::GenContext, m::DistributionMeasure, sz::Dims, ::CPUnit) =
    _flat_powrand(get_rng(ctx), get_precision(ctx), m.obj, sz)
_dist_batched_rand(ctx::GenContext, m::DistributionMeasure, sz::Dims, ::AbstractComputeUnit) =
    MeasureBase._rand_default(ctx, m, sz, MeasureBase._NoRandImpl())

# A single variate for zero batch dimensions, flat batches otherwise:
_flat_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution, ::Tuple{}) where {T<:Real} = convert_realtype(T, rand(rng, d))
_flat_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution, sz::Dims) where {T<:Real} = _flat_powrand_batch(rng, T, d, sz)

function _flat_powrand_batch(rng::AbstractRNG, ::Type{T}, d::Distribution{<:ArrayLikeVariate{0}}, sz::Dims) where {T<:Real}
    convert_realtype(T, reshape(rand(rng, d, prod(sz)), sz...))
end

function _flat_powrand_batch(rng::AbstractRNG, ::Type{T}, d::Distribution{<:ArrayLikeVariate{1}}, sz::Dims) where {T<:Real}
    convert_realtype(T, reshape(rand(rng, d, prod(sz)), size(d)..., sz...))
end

function _flat_powrand_batch(rng::AbstractRNG, ::Type{T}, d::ReshapedDistribution{N,<:Any,<:Distribution{<:ArrayLikeVariate{1}}}, sz::Dims) where {T<:Real,N}
    convert_realtype(T, reshape(rand(rng, d.dist, prod(sz)), d.dims..., sz...))
end

function _flat_powrand_batch(rng::AbstractRNG, ::Type{T}, d::Distribution, sz::Dims) where {T<:Real}
    flatview(ArrayOfSimilarArrays(convert_realtype(T, rand(rng, d, sz))))
end



@inline DensityInterface.densityof(m::DistributionMeasure) = densityof(m.obj)
@inline DensityInterface.logdensityof(m::DistributionMeasure) = logdensityof(m.obj)

@inline MeasureBase.logdensity_def(m::DistributionMeasure, x) = DensityInterface.logdensityof(m.obj, x)

# Distributions evaluate flat batches of array variates (the trailing
# dimensions are batch dimensions) directly, univariate wrappers broadcast
# their point kernels:
for (bhead, phead) in ((:batched_logdensityof_impl, :logdensityof_impl), (:batched_logdensity_def, :logdensity_def))
    @eval function MeasureBase.$bhead(m::DistributionMeasure{<:ArrayLikeVariate{N}}, X::AbstractArray) where {N}
        Distributions.logpdf(m.obj, X)
    end
    @eval function MeasureBase.$bhead(m::DistributionMeasure{<:ArrayLikeVariate{0}}, X::AbstractArray)
        MeasureBase._scalar_kernel_broadcast(MeasureBase.$phead, m, X)
    end
end
@inline MeasureBase.unsafe_logdensityof(m::DistributionMeasure, x) = DensityInterface.logdensityof(m.obj, x)
@inline MeasureBase.insupport(m::DistributionMeasure, x) = Distributions.insupport(m.obj, x) & _finite_variate(m.obj, x)
# Infinite values lie outside the support of univariate distributions,
# where Distributions may evaluate to NaN:
@inline _finite_variate(::Distribution{Univariate}, x) = isfinite(x)
@inline _finite_variate(::Distribution, x) = true

@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate{0},<:Continuous}) = Lebesgue()
@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate,<:Continuous}) = Lebesgue()^size(m.obj)
@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate{0},<:Discrete}) = Counting()
@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate,<:Discrete}) = Counting()^size(m.obj)

@inline MeasureBase.basemeasure(m::DistributionMeasure) = rootmeasure(m)

@inline MeasureBase.massof(::DistributionMeasure) = static(1.0)

@inline MeasureBase.mspace_elsize(d::Distribution) = MeasureBase.NoMSpaceElementSize{typeof(d)}()
@inline MeasureBase.mspace_flatsize(d::Distribution) = MeasureBase.NoMSpaceElementSize{typeof(d)}()
@inline MeasureBase.mspace_elsize(d::Distribution{Univariate}) = ()
@inline MeasureBase.mspace_elsize(d::Distribution{<:ArrayLikeVariate}) = size(d)
@inline MeasureBase.mspace_flatsize(d::Distribution{Univariate}) = ()
@inline MeasureBase.mspace_flatsize(d::Distribution{<:ArrayLikeVariate}) = size(d)
@inline MeasureBase.mspace_elsize(m::DistributionMeasure) = MeasureBase.mspace_elsize(m.obj)
@inline MeasureBase.mspace_flatsize(m::DistributionMeasure) = MeasureBase.mspace_flatsize(m.obj)
@inline MeasureBase.mspace_flatsize(::Type{<:Distribution{Univariate}}) = ()
@inline MeasureBase.mspace_ndims(::Type{<:Distribution{<:ArrayLikeVariate{N}}}) where {N} = N
@inline MeasureBase.mspace_ndims(::Type{AsMeasure{D}}) where {D<:Distribution} = MeasureBase.mspace_ndims(D)
@inline MeasureBase.mspace_flatsize(::Type{AsMeasure{D}}) where {D<:Distribution} = MeasureBase.mspace_flatsize(D)

@inline MeasureBase.preferred_stdmeasure(::Type{AsMeasure{D}}) where {D<:Distribution} = MeasureBase.preferred_stdmeasure(D)

@inline MeasureBase.getdof(m::DistributionMeasure{<:ArrayLikeVariate{0}}) = 1

# Delegate transport to the wrapped distribution:
@inline MeasureBase.transport_to_std(::Type{S}, m::DistributionMeasure, x) where {S<:StdMeasure} =
    MeasureBase.transport_to_std(S, m.obj, x)
@inline MeasureBase.transport_from_std(::Type{S}, m::DistributionMeasure, z) where {S<:StdMeasure} =
    MeasureBase.transport_from_std(S, m.obj, z)
@inline MeasureBase.batched_transport_to_std(::Type{S}, m::DistributionMeasure, X::AbstractArray) where {S<:StdMeasure} =
    MeasureBase.batched_transport_to_std(S, m.obj, X)
@inline MeasureBase.batched_transport_from_std(::Type{S}, m::DistributionMeasure, Z::AbstractArray) where {S<:StdMeasure} =
    MeasureBase.batched_transport_from_std(S, m.obj, Z)

@inline MeasureBase.paramnames(m::DistributionMeasure) = propertynames(m.obj)
@inline MeasureBase.params(m::DistributionMeasure) = NamedTuple{propertynames(m.obj)}(Distributions.params(m.obj))

# @inline MeasureBase.testvalue(m::DistributionMeasure) = testvalue(basemeasure(d))
