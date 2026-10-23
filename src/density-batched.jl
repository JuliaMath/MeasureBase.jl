# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Batched-first density evaluation.
#
# A flat batch of variates is an array `(variate dims..., batch dims...)`,
# zero batch dims meaning a single variate. Nested batches (arrays of
# variates with flat storage, see ArraysOfArrays) are fused at the entry
# points, the container then determines the batch dimensions. Batched
# kernels know the variate rank of their measure (see `mspace_ndims`) and
# return arrays over the batch dimensions, a number for a single variate.
# Structural measures implement their kernels once, in terms of the kernels
# of their components.

export logdensities

"""
    logdensities(μ::AbstractMeasure, X)

Compute the log-density of `μ` at each variate in the batch `X`.

`X` is an array of variates (e.g. an `ArraysOfArrays.ArrayOfSimilarArrays`,
or an array of numbers for measures with scalar variates), the flat storage
of a batch (an array of numbers with the variate dimensions leading, see
[`MeasureBase.mspace_ndims`](@ref)), or a tuple resp. `NamedTuple` of
batches for measures with tuple resp. `NamedTuple` variates. Returns an
array over the batch dimensions (a number for a single variate),
semantically equivalent to `logdensityof.(Ref(μ), X)` for arrays of
variates. Batches with flat storage are evaluated in fused operations,
compatible with GPU and traced arrays.

Measure types implement [`MeasureBase.batched_logdensityof_impl`](@ref).
"""
function logdensities end

@inline logdensities(μ::AbstractMeasure, X) = _materialize(_batched_ld(logdensityof_impl, μ, X))

"""
    logdensities(ms::AbstractArray{<:AbstractMeasure}, X)

Compute the log-densities of the measures in `ms` at the variates of the
batch `X`, paired with the measures as in broadcasting over the batch
dimensions of `X`.

`X` is a batch as for `logdensities(μ, X)` with array or numeric
variates, its batch dimensions broadcast against the dimensions of `ms`.
A single variate (no batch dimensions) is shared by all measures. This is
equivalent to `logdensityof.(ms, Xs)` for the array `Xs` of the variates
in `X`, and to `logdensityof.(ms, Ref(x))` for a single array variate
`x`, which share the implementation.

Measures with scalar variates are evaluated in one fused broadcast, over
the parameter columns of `ms` if it is a struct array (broadcasting over
a struct array of parameters returns one), also on GPUs. Arrays of powers
of the same shape broadcast their base measures over the power
dimensions. Other measures with array variates are evaluated measure by
measure.
"""
@inline logdensities(ms::AbstractArray{<:AbstractMeasure}, X) = _materialize(_measures_ld(ms, X))

"""
    MeasureBase.batched_logdensityof_impl(μ::AbstractMeasure, X)

Log-densities of `μ` (relative to its root measure) at the variates of the
flat batch `X`, an array `(variate dims..., batch dims...)`, returned as an
array over the batch dimensions. `X` may be a single variate (zero batch
dimensions), the result is a number then. Results may be lazy broadcasts,
callers materialize them where necessary. The results must be `-Inf` for
variates outside the support of `μ`.

This is the primary density extension point. The default implementation
broadcasts the point kernel [`MeasureBase.logdensityof_impl`](@ref) over
`X` for measures with scalar variates and maps it over the variate slices
of `X` (in a host loop) for measures with array variates of a declared
number of dimensions (see [`MeasureBase.mspace_ndims`](@ref)). Measure
types with array variates should implement `batched_logdensityof_impl`
directly. Structural measures evaluate their components via
`batched_logdensityof_impl` as well.
"""
function batched_logdensityof_impl end

@inline function batched_logdensityof_impl(μ::AbstractMeasure, X)
    _default_batched_kernel(logdensityof_impl, μ, X, _static_ndims(μ))
end

# The variate rank as a static integer, from the type where known:
@inline _static_ndims(μ::MU) where {MU} = _static_ndims(mspace_ndims(MU), μ)
@inline _static_ndims(n::Integer, μ) = static(n)
@inline _static_ndims(::NoMSpaceElementSize, μ) = _static_ndims_of(mspace_ndims(μ))
@inline _static_ndims_of(n::Integer) = static(n)
@inline _static_ndims_of(n::NoMSpaceElementSize) = n

"""
    MeasureBase.batched_logdensity_def(μ::AbstractMeasure, X)

Batched form of [`logdensity_def`](@ref): log-densities of `μ` relative to
`basemeasure(μ)` at the variates of the flat batch `X`, with the same
conventions and defaults as [`MeasureBase.batched_logdensityof_impl`](@ref).
"""
function batched_logdensity_def end

@inline function batched_logdensity_def(μ::AbstractMeasure, X)
    _default_batched_kernel(logdensity_def, μ, X, _static_ndims(μ))
end

@inline _default_batched_kernel(f::F, μ, X, n::Integer) where {F} = _default_batched_kernel(f, μ, X, static(n))
@inline _default_batched_kernel(f::F, μ, X, ::StaticInteger{0}) where {F} = _scalar_kernel_broadcast(f, μ, X)
@inline _default_batched_kernel(f::F, μ, X::AbstractArray, ::StaticInteger{0}) where {F} = _scalar_kernel_broadcast(f, μ, X)
@inline function _default_batched_kernel(f::F, μ, X::AbstractArray, ::StaticInteger{K}) where {F,K}
    _map_variate_slices(f, μ, X, Val(K))
end
@noinline function _default_batched_kernel(f::F, μ, X, ::NoMSpaceElementSize) where {F}
    throw(ArgumentError("Batched density evaluation for measures of type $(nameof(typeof(μ))) requires MeasureBase.mspace_ndims to be declared for the type or MeasureBase.batched_logdensityof_impl to be implemented"))
end

# Point kernels of scalar-variate measures broadcast over the batch. Static
# results are made dynamic to keep reductions type stable.
@inline function _scalar_kernel_broadcast(f::F, μ, X::AbstractArray) where {F}
    Broadcast.instantiate(Broadcast.broadcasted(_DynamicLogd(f, μ), X))
end
@inline _scalar_kernel_broadcast(f::F, μ, x) where {F} = _dynamic_logd(f(μ, x), x)

struct _DynamicLogd{F,M} <: Function
    f::F
    μ::M
end
@inline (k::_DynamicLogd)(x) = _dynamic_logd(k.f(k.μ, x), x)

# Point kernels of array-variate measures map over the variate slices of
# the batch, a batch of a single variate is evaluated directly:
@inline _map_variate_slices(f::F, μ, X::AbstractArray{<:Any,K}, ::Val{K}) where {F,K} = _dynamic_logd(f(μ, X), X)
@inline function _map_variate_slices(f::F, μ, X::AbstractArray, ::Val{K}) where {F,K}
    map(_DynamicLogd(f, μ), sliced(X, Val(K)))
end

# The batched kernel for a point-level density function:
@inline _batched_kernel(::typeof(logdensityof_impl), μ, X) = batched_logdensityof_impl(μ, X)
@inline _batched_kernel(::typeof(logdensity_def), μ, X) = batched_logdensity_def(μ, X)

# Flat storage of a (nested) batch: the underlying array of memory-ordered
# split arrays, a stacked copy for other known split modes. Struct arrays
# of tuple and named tuple variates have the flat storage of their
# components (copied where necessary, the batch dimensions are shared).
struct NoFlatStorage end
@inline _flat_storage(X::AbstractArray{<:Number}) = X
@inline _flat_storage(X::AbstractArray) = _flat_storage_bymode(X, getsplitmode(X))
# Static arrays of static arrays stack into static arrays (their
# reinterpreted flat views have inconsistent axes):
@inline function _flat_storage(X::StaticArray{<:Tuple,<:StaticArray{<:Tuple,<:Number}})
    SArray{Tuple{size(eltype(X))...,size(X)...}}(mapreduce(Tuple, (a, b) -> (a..., b...), X))
end
@inline _flat_storage(x) = NoFlatStorage()
@inline function _flat_storage_bymode(X::AbstractArray, smode::AbstractSplitMode)
    _flat_storage(is_memordered_splitmode(smode) ? fused(X) : stacked(X))
end
@inline _flat_storage_bymode(::AbstractArray, ::UnknownSplitMode) = NoFlatStorage()
@inline _flat_storage_bymode(::AbstractArray, ::NonSplitMode) = NoFlatStorage()
# Ragged batches (e.g. `VectorOfArrays`) are evaluated variate by variate:
@inline _flat_storage_bymode(::AbstractArray, ::AbstractPartMode) = NoFlatStorage()

@inline _flat_storage(X::StructArray{<:Union{Tuple,NamedTuple}}) = _components_storage(StructArrays.components(X))
@inline function _flat_storage(X::AbstractArray{<:Union{Tuple,NamedTuple}})
    _flat_storage(StructArray(X; unwrap = T -> T <: Union{Tuple,NamedTuple}))
end
@inline _components_storage(cs::Tuple) = map(_component_storage, cs)
@inline _components_storage(cs::NamedTuple{names}) where {names} = NamedTuple{names}(map(_component_storage, values(cs)))
@inline _component_storage(c::StructArray{<:Union{Tuple,NamedTuple}}) = _flat_storage(c)
@inline _component_storage(c::AbstractArray{<:Number}) = c
@inline _component_storage(c::AbstractArray{<:AbstractArray}) = _component_flat(c, _flat_storage(c))
@inline _component_flat(c, c_flat::AbstractArray) = c_flat
@inline _component_flat(c, ::NoFlatStorage) = stacked(c)

# Entry: arrays of numbers are flat storage, arrays of variates are fused
# into their flat storage (else evaluated variate by variate), tuples and
# named tuples of batches go to the kernels directly.
@inline _batched_ld(f::F, μ, X::AbstractArray{<:Number}) where {F} = _checked_batched_kernel(f, μ, X)
@inline _batched_ld(f::F, μ, x::Number) where {F} = _point_ld(f, μ, x)
@inline _batched_ld(f::F, μ, X::Union{Tuple,NamedTuple}) where {F} = _batched_kernel(f, μ, X)
@inline _batched_ld(f::F, μ, X::AbstractArray) where {F} = _batched_ld_nested(f, μ, X, _flat_storage(X), _static_ndims(μ))
@inline function _batched_ld_nested(f::F, μ, X::AbstractArray, X_flat::AbstractArray, ::StaticInteger) where {F}
    _check_batch_shape(_checked_batched_kernel(f, μ, X_flat), X)
end
@inline function _batched_ld_nested(f::F, μ, X::AbstractArray, X_flat::Union{Tuple,NamedTuple}, ::Any) where {F}
    _check_batch_shape(_batched_kernel(f, μ, X_flat), X)
end
@inline function _batched_ld_nested(f::F, μ, X::AbstractArray, ::Any, ::Any) where {F}
    Broadcast.instantiate(Broadcast.broadcasted(_PointLogd(f, μ), X))
end

struct _PointLogd{F,M} <: Function
    f::F
    μ::M
end
@inline (k::_PointLogd)(x) = _point_ld(k.f, k.μ, x)

# Kernels rely on the leading dimensions of flat batches to match the flat
# variate size of the measure where it is known, singleton dimensions would
# broadcast. The entry points check it, and so do structural measures
# without a flat size of their own when they pass batches to components:
@inline function _checked_batched_kernel(f::F, μ, X::AbstractArray{<:Number}) where {F}
    _check_flatsize(X, mspace_flatsize(μ))
    _batched_kernel(f, μ, X)
end
@inline _checked_batched_kernel(f::F, μ, X) where {F} = _batched_kernel(f, μ, X)

# Results over an array of variates must have the shape of the array, a
# mismatch means the variate rank of the measure doesn't match the batch:
@inline function _check_batch_shape(result, X::AbstractArray)
    if size(result) != size(X)
        _throw_batch_shape(size(result), size(X))
    end
    return result
end
@noinline function _throw_batch_shape(sz_result, sz_batch)
    throw(ArgumentError("Batched density kernel returned a result of size $sz_result for a batch of size $sz_batch, the variate dimensions of the measure don't match the batch"))
end

# Point evaluation: array variates of measures with a declared variate rank
# go through the batched kernel with zero batch dimensions, other variates
# through the point kernel. A batched kernel that returns an array for a
# single variate has taken the variate for a batch: the variate doesn't
# fit the measure, or the measure lacks a batched kernel for array
# variates.
@inline _point_ld(f::F, μ, x::AbstractArray{<:Number}) where {F} = _point_ld_byrank(f, μ, x, _static_ndims(μ))
@inline _point_ld(f::F, μ, x) where {F} = f(μ, x)
@inline _point_ld_byrank(f::F, μ, x, ::StaticInteger) where {F} = _point_result(_materialize(_checked_batched_kernel(f, μ, x)), μ)
@inline _point_ld_byrank(f::F, μ, x, ::NoMSpaceElementSize) where {F} = f(μ, x)
@inline _point_ld(f::F, μ::PrimitiveMeasure, x::AbstractArray{<:Number}) where {F} = f(μ, x)

@inline _point_result(ℓ::Number, μ) = ℓ
@inline _point_result(ℓ::AbstractArray{<:Number,0}, μ) = ℓ[]
@noinline function _point_result(ℓ, μ)
    throw(ArgumentError("Density evaluation of measures of type $(nameof(typeof(μ))) at an array variate resulted in a batch of densities: the variate doesn't fit the measure, or the measure lacks a batched kernel for array variates"))
end

const _LazyBroadcast = Broadcast.Broadcasted

@noinline _throw_size_mismatch() = throw(ArgumentError("Size of variate doesn't match size of measure"))

# The leading dimensions of a flat batch must match a flat variate size:
@inline function _check_flatsize(A::AbstractArray, sz_flat::SizeLike)
    n = length(size_tuple(sz_flat))
    if ndims(A) < n || ntuple(i -> size(A, i), Val(n)) != asnonstatic(sz_flat)
        _throw_size_mismatch()
    end
    return nothing
end
@inline _check_flatsize(::AbstractArray, ::NoMSpaceElementSize) = nothing

@inline _materialize(bc::_LazyBroadcast) = copy(bc)
@inline _materialize(x) = x

# `sum_leading_dims` of log-densities. Partial sums of lazy broadcasts over
# host arrays don't materialize the elementwise log-densities (Base can't
# reduce selected dimensions of a `Broadcasted`):
@inline _sum_leading_logd(ℓ, n::StaticInteger; kw...) = sum_leading_dims(ℓ, n; kw...)
function _sum_leading_logd(ℓ::Broadcast.Broadcasted{<:Broadcast.DefaultArrayStyle}, n::StaticInteger{N}; kw...) where {N}
    axs = axes(ℓ)
    (0 < N < length(axs) && !isempty(ℓ)) || return sum_leading_dims(ℓ, n; kw...)
    lead = CartesianIndices(ntuple(i -> axs[i], Val(N)))
    map(CartesianIndices(ntuple(i -> axs[N + i], Val(length(axs) - N)))) do j
        s = zero(ℓ[CartesianIndex(first(lead), j)])
        @inbounds @simd for i in lead
            s += ℓ[CartesianIndex(i, j)]
        end
        s
    end
end


@inline _lazy_add(a::Number, b::Number) = a + b
@inline _lazy_add(a, b) = Broadcast.instantiate(Broadcast.broadcasted(+, a, b))

# Zero log-densities over the batch dimensions of a flat batch of variates
# with `n` variate dimensions:
@inline function _zero_logd_batch(X::AbstractArray, n::Integer)
    FillArrays.Zeros{_logd_numtype(X)}(ntuple(i -> size(X, n + i), ndims(X) - n))
end
@inline function _zero_logd_batch(X::AbstractArray{<:Any,N}, ::StaticInteger{K}) where {N,K}
    FillArrays.Zeros{_logd_numtype(X)}(ntuple(i -> size(X, K + i), Val(N - K)))
end
@inline _zero_logd_batch(X::AbstractArray{<:Any,N}, ::StaticInteger{N}) where {N} = zero(_logd_numtype(X))
@inline _zero_logd_batch(x::Number, ::StaticInteger{0}) = zero(_logd_numtype(x))
@noinline function _zero_logd_batch(X, ::NoMSpaceElementSize{MU}) where {MU}
    throw(ArgumentError("Batched density evaluation for measures of type $(nameof(MU)) requires MeasureBase.mspace_ndims to be declared for the type or MeasureBase.batched_logdensity_def to be implemented"))
end


# Streams: variates of composed measures are consumed from flat vector
# streams, batched as `(rows, batch dims...)`.

"""
    MeasureBase.batched_logdensityof_with_rest(μ::AbstractMeasure, X, sz::SizeLike)

Consume variates of `μ` from the batch `X` of flat vector streams (first
dimension along the streams, further dimensions are batch dimensions), a
batch of variates of size `sz` per stream, and compute their
log-densities.

Returns a tuple `(ℓ, X_rest)` of the log-densities, an array over
`(sz..., batch dims...)` (possibly lazy, a number for a single stream and
`sz == ()`), and the unconsumed rest of the streams. Measure types whose
variates have a fixed size consume `prod(sz)` variates in one batched
kernel evaluation, the default implementation does so for the variate size
given by [`MeasureBase.mspace_flatsize`](@ref) or
[`MeasureBase.some_mspace_elsize`](@ref). Measure types with variates of
value-dependent size implement `batched_logdensityof_with_rest` for
single streams and `sz == ()` themselves and report
`MeasureBase.fixed_stream_size` as false, so that the enclosing stream
combinators consume batches stream by stream.
"""
function batched_logdensityof_with_rest end

function batched_logdensityof_with_rest(μ::AbstractMeasure, X::AbstractArray, sz::SizeLike)
    _stream_ld_with_rest(logdensityof_impl, μ, X, sz)
end

# A single stream consumes one variate via the point path:
function batched_logdensityof_with_rest(μ::AbstractMeasure, x::AbstractVector, ::Tuple{})
    ℓ, _, x_rest = logdensityof_with_rest(μ, x)
    return ℓ, x_rest
end

function _stream_ld_with_rest(f::F, μ, X::AbstractArray, sz::SizeLike) where {F}
    vsz = _stream_consume_size(μ)
    X_μ, X_rest = _batched_consume(X, vsz, sz)
    return _consumed_ld(f, μ, X_μ, vsz), X_rest
end

# Scalar variates are consumed as `(1, sz..., batch dims...)`, their
# leading dimension is dropped before the kernel runs:
@inline _consumed_ld(f::F, μ, X_μ, ::Tuple{}) where {F} = _batched_kernel(f, μ, _drop_stdstream_dim(X_μ))
@inline _consumed_ld(f::F, μ, X_μ, ::SizeLike) where {F} = _batched_kernel(f, μ, X_μ)

# Consume `prod(sz)` variates of flat size `vsz` from the leading rows of a
# batch of streams as a flat batch `(vsz..., sz..., batch dims...)`; scalar
# variates as `(1, sz..., batch dims...)`. Static sizes keep static
# streams static.
@inline function _batched_consume(X::AbstractArray, vsz::SizeLike, sz::SizeLike)
    dims = _consumed_dims(vsz)
    X_flat, X_rest = _batched_split(X, _chunk_rows(prod(dims), sz))
    return _reshape_consumed(X_flat, (dims..., size_tuple(sz)...)), X_rest
end
@inline _consumed_dims(::Tuple{}) = (static(1),)
@inline _consumed_dims(vsz::SizeLike) = size_tuple(vsz)
@inline _chunk_rows(n::IntegerLike, sz::SizeLike) = n * size2length(sz)

@inline _reshape_consumed(X_flat::AbstractArray, ::Tuple{IntegerLike}) = X_flat
@inline function _reshape_consumed(X_flat::AbstractArray, dims::Tuple{Vararg{IntegerLike}})
    maybestatic_reshape(X_flat, (dims..., Base.tail(_batch_dims(X_flat))...))
end

# Sizes as tuples of (maybe static) integers and the leading dimension of
# a batch of streams:
@inline _batch_dims(A::AbstractArray) = size_tuple(maybestatic_size(A))
@inline _as_stdstream_batch(Z::AbstractArray) = merge_leading_dims(Z, static(0))
@inline _as_stdstream_batch(z::Number) = SVector(z)
@inline _drop_stdstream_dim(Z::AbstractArray) = drop_leading_dims(Z, static(1))

@inline function _batched_split(A::AbstractArray, n::IntegerLike)
    n_rows = dynamic(n)
    stream_idxs = axes(A, 1)
    if length(stream_idxs) < n_rows
        throw(ArgumentError("Variate streams too short during batched evaluation"))
    end
    batch_axes = Base.tail(axes(A))
    i0 = first(stream_idxs)
    A_flat = view(A, i0:(i0 + n_rows - 1), batch_axes...)
    A_rest = view(A, (i0 + n_rows):last(stream_idxs), batch_axes...)
    return A_flat, A_rest
end

# Static streams split into static chunks for static row counts:
@inline _batched_split(A::StaticArray, n_rows::StaticInteger) = maybestatic_split(A, n_rows)

@noinline function _throw_stream_too_long()
    throw(ArgumentError("Variate streams too long during density evaluation"))
end

"""
    MeasureBase.fixed_stream_size(::Type{MU})::StaticBool
    MeasureBase.fixed_stream_size(μ)::StaticBool

Whether the variates that measures of type `MU` consume from streams have
the same size for all streams, so that batches of streams can be consumed
in fused operations. Otherwise the outermost stream combinator consumes a
batch of streams stream by stream.

Defaults to whether the variate rank of `MU` is known (see
[`MeasureBase.mspace_ndims`](@ref)). Measure types whose variate sizes
depend on variate values (e.g. binds) must return `static(false)`.
"""
function fixed_stream_size end

@inline fixed_stream_size(μ::MU) where {MU} = fixed_stream_size(MU)
@inline fixed_stream_size(::Type{MU}) where {MU} = static(mspace_ndims(MU) isa Integer)

# Batches of streams consumed stream by stream (host loop):
function _streamwise_ld(::typeof(logdensityof_impl), μ, X::AbstractArray)
    map(sliced(X, Val(1))) do x
        ℓ, x_rest = batched_logdensityof_with_rest(μ, x, ())
        isempty(x_rest) || _throw_stream_too_long()
        _dynamic_logd(_materialize(ℓ), x)
    end
end


# Parameter batches: arrays of measures at batches of variates, paired as
# in broadcasting.

@inline _measures_ld(ms::AbstractArray{M}, X) where {M} = _measures_ld(ms, X, _static_ndims_of(mspace_ndims(M)))

# Scalar variates in one broadcast over the measures (their parameter
# columns for struct arrays) and the batch:
@inline _measures_ld(ms, x::Number, ::Any) = _measures_bc(ms, x)
@inline _measures_ld(ms, X::AbstractArray{<:Number}, ::StaticInteger{0}) = _measures_bc(ms, X)
@inline function _measures_ld(ms, X::AbstractArray{<:Number}, k::StaticInteger)
    ndims(X) >= k || _throw_size_mismatch()
    _measures_array_ld(ms, X, k)
end
# Array variates measure by measure, a single variate is shared by all:
@inline _measures_array_ld(ms, X::AbstractArray, ::StaticInteger{K}) where {K} = _measures_bc(ms, _variates_arg(X, Val(K)))
@inline _variates_arg(X::AbstractArray{<:Any,K}, ::Val{K}) where {K} = Ref(X)
@inline _variates_arg(X::AbstractArray, ::Val{K}) where {K} = sliced(X, Val(K))
# Arrays of variates via their flat storage where they have one:
@inline _measures_ld(ms, X::AbstractArray, k) = _measures_ld_nested(ms, X, _flat_storage(X), k)
@inline _measures_ld_nested(ms, ::AbstractArray, X_flat::AbstractArray{<:Number}, k::StaticInteger) = _measures_ld(ms, X_flat, k)
@inline _measures_ld_nested(ms, X::AbstractArray, ::Any, ::Any) = _measures_bc(ms, X)
@noinline function _measures_ld(ms, X, ::Any)
    throw(ArgumentError("logdensities over arrays of measures supports numeric and array variates, broadcast logdensityof over arrays of other variates"))
end

@inline _measures_bc(ms, X) = _marginal_broadcast(_DynamicPointLogd(logdensityof), ms, X)

# Broadcasting log-densities over arrays of measures shares the
# implementation of `logdensities`:
@inline function Broadcast.broadcasted(::typeof(logdensityof), ms::AbstractArray{<:AbstractMeasure}, X)
    _measures_bc(ms, X)
end
@inline function Broadcast.broadcasted(::typeof(logdensityof), ms::AbstractArray{M}, x::Base.RefValue{<:AbstractArray{<:Number}}) where {M<:AbstractMeasure}
    _shared_variate_ld(ms, x, _static_ndims_of(mspace_ndims(M)))
end
@inline function _shared_variate_ld(ms, x::Base.RefValue{<:AbstractArray{<:Number,K}}, k::StaticInteger{K}) where {K}
    _measures_ld(ms, x[], k)
end
@inline _shared_variate_ld(ms, x::Base.RefValue, ::Any) = _measures_bc(ms, x)

@inline function Broadcast.broadcasted(::typeof(densityof), ms::AbstractArray{<:AbstractMeasure}, X)
    Broadcast.broadcasted(exp, Broadcast.broadcasted(logdensityof, ms, X))
end
