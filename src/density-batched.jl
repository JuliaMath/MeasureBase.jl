# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

export logdensities

"""
    logdensities(μ::AbstractMeasure, X::AbstractArray)

Compute the log-density of `μ` at each point in `X`.

Returns an array of the shape of `X`, semantically equivalent to
`logdensityof.(Ref(μ), X)`. Batches with flat storage of a known variate
size (e.g. `ArraysOfArrays.ArrayOfSimilarArrays`) are evaluated in one
fused operation over the flat data, compatible with GPU-backed storage.

For measures with array-valued variates, `X` may also be the flat storage
of the batch itself, with the variate dimensions leading (see
[`MeasureBase.mspace_flatsize`](@ref)). The result then has the remaining
dimensions of `X`.

Measure types should specialize
[`MeasureBase.batched_logdensityof_impl`](@ref) instead of `logdensities`
itself.
"""
function logdensities end

@inline logdensities(μ::AbstractMeasure, X::AbstractArray) = _materialize(_batched_ld(logdensityof_impl, μ, X))

"""
    MeasureBase.batched_logdensityof_impl(μ::AbstractMeasure, A::AbstractArray)

Implements [`logdensities`](@ref) for a batch `A` of variates of `μ` in
flat storage: the leading dimensions of `A` are the variate dimensions
(see [`MeasureBase.mspace_flatsize`](@ref)), any further dimensions are
batch dimensions. Returns the log-densities as an array over the batch
dimensions, or a scalar if there are none.

Power measures never reach `batched_logdensityof_impl`, their power
structure is unwrapped beforehand. Implementations must handle points
outside the support of `μ` (the result must be `-Inf` there).

The default implementation broadcasts the log-density over `A` for
measures with scalar variates and maps it over the variate slices of `A`
otherwise. The result may be a lazy broadcast, callers materialize it
where necessary.
"""
function batched_logdensityof_impl end

@inline function batched_logdensityof_impl(μ::AbstractMeasure, A::AbstractArray)
    _batched_ld_generic(logdensityof_impl, μ, A)
end

# Batched density machinery, parameterized over the point-level density
# function `f` (`logdensityof_impl` or `logdensity_def`).
#
# Variates and batches with flat storage are evaluated in one call of the
# batched point kernel of the base measure: the leading dimensions of the
# flat data are the variate dimensions of the base measure, followed by the
# power dimensions and any batch dimensions. The power dimensions are summed
# afterwards. Nested arrays without flat storage of a known layout are
# evaluated level by level.

struct NoFlatStorage end

@inline _batched_ld(f::F, μ, X::AbstractArray) where {F} = _batched_ld_sized(f, μ, X, mspace_flatsize(μ))

@inline function _batched_ld_sized(f::F, μ, X::AbstractArray, sz_flat::SizeLike) where {F}
    _batched_ld_flat(f, μ, X, _flat_storage(X), sz_flat)
end

@inline function _batched_ld_sized(f::F, μ, X::AbstractArray, ::NoMSpaceElementSize) where {F}
    map(x -> _pointwise_ld(f, μ, x), X)
end

@inline function _batched_ld_flat(f::F, μ, X, X_flat::AbstractArray, sz_flat) where {F}
    ν, n_pwr = _pwr_unwrap(μ)
    _check_flatsize(X_flat, sz_flat)
    _sum_leading_dims(_batched_kernel(f, ν, X_flat), n_pwr)
end

@inline function _batched_ld_flat(f::F, μ, X, ::NoFlatStorage, sz_flat) where {F}
    map(x -> _pointwise_ld(f, μ, x), X)
end

@inline _pointwise_ld(f::F, μ, x) where {F} = f(μ, x)
@inline _pointwise_ld(f::F, μ::PowerMeasure, x) where {F} = _powered_ld(f, μ, x)

# Log-density of a power measure at a single variate:

@inline _powered_ld(f::F, μ::PowerMeasure, x) where {F} = _powered_ld_sized(f, μ, x, mspace_flatsize(μ))

@inline function _powered_ld_sized(f::F, μ::PowerMeasure, x, sz_flat::SizeLike) where {F}
    _powered_ld_flat(f, μ, x, _flat_storage(x), sz_flat)
end

@inline function _powered_ld_sized(f::F, μ::PowerMeasure, x, ::NoMSpaceElementSize) where {F}
    _powered_ld_pointwise(f, μ, x)
end

@inline function _powered_ld_flat(f::F, μ::PowerMeasure, x, x_flat::AbstractArray, sz_flat) where {F}
    ν, n_pwr = _pwr_unwrap(μ)
    _check_flatsize(x_flat, sz_flat)
    _sum_leading_dims(_batched_kernel(f, ν, x_flat), n_pwr)
end

@inline _powered_ld_flat(f::F, μ::PowerMeasure, x, ::NoFlatStorage, sz_flat) where {F} = _powered_ld_pointwise(f, μ, x)

# Sum of the point-level densities over the elements of the variate:
@inline function _powered_ld_pointwise(f::F, μ::PowerMeasure, x::AbstractArray) where {F}
    if maybestatic_size(x) != pwr_size(μ)
        throw(ArgumentError("Size of variate doesn't match size of power measure"))
    end
    sum(Base.Fix1(_pointwise_ld_dyn, (f, pwr_base(μ))), x)
end

function _powered_ld_pointwise(f::F, ::PowerMeasure, x) where {F}
    throw(ArgumentError("Variate of a power measure must be an array"))
end

@inline _pointwise_ld_dyn((f, μ), x) = dynamic(_pointwise_ld(f, μ, x))

@inline _pwr_unwrap(μ) = (μ, static(0))
@inline function _pwr_unwrap(μ::PowerMeasure)
    ν, n = _pwr_unwrap(pwr_base(μ))
    ν, n + static(length(pwr_axes(μ)))
end

# Flat storage of a (nested) variate or batch: the underlying array of
# memory-ordered split arrays, a stacked copy for other known split modes.
# Nested arrays of unknown layout have no flat storage.

@inline _flat_storage(x::AbstractArray{<:Number}) = x
@inline _flat_storage(x::AbstractArray) = _flat_storage_bymode(x, getsplitmode(x))
@inline _flat_storage(x) = NoFlatStorage()

@inline function _flat_storage_bymode(x::AbstractArray, smode::AbstractSplitMode)
    _flat_storage(is_memordered_splitmode(smode) ? fused(x) : stacked(x))
end
@inline _flat_storage_bymode(::AbstractArray, ::UnknownSplitMode) = NoFlatStorage()
@inline _flat_storage_bymode(::AbstractArray, ::NonSplitMode) = NoFlatStorage()

@inline function _check_flatsize(A::AbstractArray, sz_flat::SizeLike)
    n = length(sz_flat)
    if ndims(A) < n || ntuple(i -> size(A, i), Val(n)) != Tuple(sz_flat)
        throw(ArgumentError("Size of variate doesn't match size of measure"))
    end
    return nothing
end

# Point kernel of the base measure over the flat batch:

@inline _batched_kernel(::typeof(logdensityof_impl), ν, A::AbstractArray) = batched_logdensityof_impl(ν, A)
@inline _batched_kernel(f::F, ν, A::AbstractArray) where {F} = _batched_ld_generic(f, ν, A)

# Powers of primitive measures have log-density zero relative to their base:
@inline function _batched_kernel(::typeof(logdensity_def), ν::PrimitiveMeasure, A::AbstractArray)
    FillArrays.Zeros{Float64}(size(A))
end

@inline _batched_ld_generic(f::F, ν, A::AbstractArray) where {F} = _batched_ld_byflatsize(f, ν, A, mspace_flatsize(ν))

# Scalar variates: one lazy broadcast over the whole batch, so that
# reductions over it don't need to allocate the intermediate result.
# Static results of point kernels are made dynamic, to keep reductions
# over them type stable.
@inline function _batched_ld_byflatsize(f::F, ν, A::AbstractArray, ::Tuple{}) where {F}
    Broadcast.instantiate(Broadcast.broadcasted(dynamic ∘ Base.Fix1(f, ν), A))
end

# Array variates: map over the variate slices, or evaluate directly if `A`
# is a single variate.
@inline function _batched_ld_byflatsize(f::F, ν, A::AbstractArray, sz::SizeLike) where {F}
    _batched_ld_slices(f, ν, A, Val(length(sz)))
end

# Variates of unknown size: the elements of `A` are the variates.
@inline function _batched_ld_byflatsize(f::F, ν, A::AbstractArray, ::NoMSpaceElementSize) where {F}
    map(Base.Fix1(f, ν), A)
end

@inline _batched_ld_slices(f::F, ν, A::AbstractArray{<:Any,N}, ::Val{N}) where {F,N} = f(ν, A)
@inline function _batched_ld_slices(f::F, ν, A::AbstractArray, ::Val{M}) where {F,M}
    map(Base.Fix1(f, ν), sliced(A, Val(M)))
end

# Sum over the leading `N` dimensions; a full reduction yields a scalar.
# Lazy broadcasts are reduced without materialization where the broadcast
# style supports it, and materialized before partial reductions.

const _LazyBroadcast = Broadcast.Broadcasted
const _EagerReducibleBroadcast = Broadcast.Broadcasted{<:Union{Broadcast.DefaultArrayStyle,StaticArrays.StaticArrayStyle}}

@inline _materialize(bc::_LazyBroadcast) = copy(bc)
@inline _materialize(x) = x

@inline _sum_leading_dims(x::Number, ::StaticInteger{0}) = x
@inline _sum_leading_dims(A::AbstractArray, n::StaticInteger) = _sum_leading_dims_impl(A, n, static(ndims(A)))
@inline _sum_leading_dims(bc::_LazyBroadcast, n::StaticInteger) = _sum_leading_dims_lazy(bc, n, static(ndims(bc)))
@inline _sum_leading_dims_lazy(bc::_LazyBroadcast, ::StaticInteger{0}, ::StaticInteger) = bc
@inline _sum_leading_dims_lazy(bc::_LazyBroadcast, ::StaticInteger{0}, ::StaticInteger{0}) = bc
@inline _sum_leading_dims_lazy(bc::_EagerReducibleBroadcast, ::StaticInteger{0}, ::StaticInteger{0}) = bc
@inline function _sum_leading_dims_lazy(bc::_EagerReducibleBroadcast, ::StaticInteger{N}, ::StaticInteger{N}) where {N}
    # Empty broadcasts have no known element type to reduce over lazily:
    isempty(bc) ? sum(copy(bc)) : sum(bc)
end
@inline _sum_leading_dims_lazy(bc::_LazyBroadcast, ::StaticInteger{N}, ::StaticInteger{N}) where {N} = sum(copy(bc))
@inline function _sum_leading_dims_lazy(bc::_LazyBroadcast, n::StaticInteger, ::StaticInteger)
    _sum_leading_dims(copy(bc), n)
end
@inline _sum_leading_dims_impl(A::AbstractArray, ::StaticInteger{0}, ::StaticInteger) = A
@inline _sum_leading_dims_impl(A::AbstractArray, ::StaticInteger{0}, ::StaticInteger{0}) = A
@inline _sum_leading_dims_impl(A::AbstractArray, ::StaticInteger{N}, ::StaticInteger{N}) where {N} = sum(A)
@inline function _sum_leading_dims_impl(A::AbstractArray, ::StaticInteger{N}, ::StaticInteger) where {N}
    dropdims(_sum_dims_seq(A, static(N)); dims = ntuple(identity, Val(N)))
end
@inline _sum_dims_seq(A::AbstractArray, ::StaticInteger{0}) = A
@inline function _sum_dims_seq(A::AbstractArray, ::StaticInteger{N}) where {N}
    _sum_dims_seq(sum(A; dims = N), static(N - 1))
end
