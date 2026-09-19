# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Batched-first transport over flat batches of variates `(variate dims...,
# batch dims...)`, zero batch dims meaning a single variate. Streams of
# standard variates are batches `(dof, batch dims...)`, consumed along their
# first dimension. Kernels know the variate rank of their measure (see
# `mspace_ndims`), structural measures implement them once in terms of the
# kernels of their components.

"""
    MeasureBase.batched_transport_to_std(::Type{S}, μ, X)

Transport the flat batch `X` of variates of `μ` to a batch `(getdof(μ),
batch dims...)` of variates of the standard measure type `S`. `X` may be a
single variate, the result is a vector then.

The default implementation broadcasts the point transport
[`MeasureBase.transport_to_std`](@ref) for measures with scalar variates
and maps it over the variate slices of `X` (in a host loop) for measures
with array variates of a declared number of dimensions (see
[`MeasureBase.mspace_ndims`](@ref)). Measure types with array variates
should implement `batched_transport_to_std` directly.
"""
function batched_transport_to_std end

@inline function batched_transport_to_std(::Type{S}, μ, X) where {S<:StdMeasure}
    _batched_to_std(S, μ, X, _static_ndims(μ))
end

@inline function _batched_to_std(::Type{S}, μ, X, ::StaticInteger{0}) where {S}
    _as_stdstream_batch(broadcast(Base.Fix1(_ToStd{S}(), μ), X))
end
@inline function _batched_to_std(::Type{S}, μ, X::AbstractArray, ::StaticInteger{0}) where {S}
    _as_stdstream_batch(broadcast(Base.Fix1(_ToStd{S}(), μ), X))
end
@inline function _batched_to_std(::Type{S}, μ, X::AbstractArray, ::StaticInteger{K}) where {S,K}
    _to_std_slices(S, μ, X, Val(K))
end
@noinline function _batched_to_std(::Type{S}, μ, X, ::NoMSpaceElementSize) where {S}
    throw(ArgumentError("Batched transport requires MeasureBase.mspace_ndims to be declared for measures of type $(nameof(typeof(μ))) or MeasureBase.batched_transport_to_std to be implemented"))
end

@inline _to_std_slices(::Type{S}, μ, X::AbstractArray{<:Any,K}, ::Val{K}) where {S,K} = _as_stdstream(transport_to_std(S, μ, X))
@inline function _to_std_slices(::Type{S}, μ, X::AbstractArray, ::Val{K}) where {S,K}
    stacked(map(Base.Fix1(_ToStd{S}(), μ), sliced(X, Val(K))))
end

# Merge the leading `N` dimensions of an array into one, `N == 0` adds a
# leading dimension of size one:
@inline function _merge_leading_dims(A::AbstractArray, ::StaticInteger{N}) where {N}
    ndims(A) >= N || _throw_size_mismatch()
    dims = _batch_dims(A)
    lead = ntuple(i -> dims[i], Val(N))
    _reshape_batch(A, (prod(lead), ntuple(i -> dims[N + i], Val(length(dims) - N))...))
end
@inline _merge_leading_dims(A::AbstractArray, ::StaticInteger{0}) = _reshape_batch(A, (static(1), _batch_dims(A)...))

# A flat batch of variates of `μ` as a batch of streams, the variate
# dimensions merged into the first dimension. Tuples of batches (tuple
# products and their powers) interleave the rows of their components
# variate by variate.
@inline _as_stream_batch(X, μ) = _as_stream_batch(X, _static_ndims(μ))
@inline _as_stream_batch(X::AbstractArray, ::StaticInteger{K}) where {K} = _merge_leading_dims(X, static(K))
@inline _as_stream_batch(x::Number, ::StaticInteger{0}) = SVector(x)
@noinline function _as_stream_batch(X, ::NoMSpaceElementSize)
    throw(ArgumentError("Concatenating batches of variates requires MeasureBase.mspace_ndims to be declared for the measures involved"))
end
@inline function _as_stream_batch(X::Union{Tuple,NamedTuple}, μ::ProductMeasure)
    vcat(map(_as_stream_batch, values(X), values(marginals(μ)))...)
end
function _as_stream_batch(X::Union{Tuple,NamedTuple}, μ::PowerMeasure)
    ν, _ = _pwr_unwrap(μ)
    n_pwr = length(_pwr_dims(μ))
    n = prod(map(dynamic, _pwr_dims(μ)))
    parts = map(values(X), values(marginals(ν))) do Xi, m
        A = _as_stream_batch(Xi, m)
        reshape(A, (size(A, 1), n, ntuple(i -> size(A, 1 + n_pwr + i), Val(ndims(A) - 1 - n_pwr))...))
    end
    _merge_leading_dims(vcat(parts...), static(2))
end

# The standard variates of a single variate must form a vector:
@inline _single_std(z::AbstractVector) = z
@noinline function _single_std(z)
    throw(ArgumentError("Transport of a single variate resulted in a batch of standard variates, the variate doesn't fit the measure"))
end

@noinline function _throw_std_length_mismatch()
    throw(ArgumentError("Length of standard variates doesn't match the degrees of freedom of the measure"))
end


"""
    MeasureBase.batched_transport_from_std(::Type{S}, μ, Z::AbstractArray)

Transport the batch `Z` of variates of the standard measure type `S`, of
size `(getdof(μ), batch dims...)`, to a flat batch of variates of `μ`. A
single stream `Z` yields a single variate.

The default implementation broadcasts the point transport
[`MeasureBase.transport_from_std`](@ref) for measures with scalar
variates and maps it over the columns of `Z` (in a host loop) for measures
with array variates of a declared number of dimensions.
"""
function batched_transport_from_std end

@inline function batched_transport_from_std(::Type{S}, μ, Z::AbstractArray) where {S<:StdMeasure}
    _batched_from_std(S, μ, Z, _static_ndims(μ))
end

@inline _batched_from_std(::Type{S}, μ, Z::AbstractArray, ::StaticInteger{0}) where {S} = _from_std_scalar(S, μ, Z)
@inline _batched_from_std(::Type{S}, μ, Z::AbstractArray, ::StaticInteger{K}) where {S,K} = _from_std_columns(S, μ, Z)
@noinline function _batched_from_std(::Type{S}, μ, ::AbstractArray, ::NoMSpaceElementSize) where {S}
    throw(ArgumentError("Batched transport requires MeasureBase.mspace_ndims to be declared for measures of type $(nameof(typeof(μ))) or MeasureBase.batched_transport_from_std to be implemented"))
end

@inline _from_std_scalar(::Type{S}, μ, z::AbstractVector) where {S} = transport_from_std(S, μ, z[begin])
@inline function _from_std_scalar(::Type{S}, μ, Z::AbstractArray) where {S}
    broadcast(Base.Fix1(_FromStd{S}(), μ), _drop_stdstream_dim(Z))
end
@inline _from_std_columns(::Type{S}, μ, z::AbstractVector) where {S} = transport_from_std(S, μ, z)
@inline function _from_std_columns(::Type{S}, μ, Z::AbstractArray) where {S}
    stacked(map(Base.Fix1(_FromStd{S}(), μ), sliced(Z, Val(1))))
end


"""
    MeasureBase.batched_transport_to_std_with_rest(::Type{S}, μ, X::AbstractArray, sz::Dims)

Consume variates of `μ` from the batch `X` of flat vector streams (first
dimension along the streams, further dimensions are batch dimensions), a
batch of variates of size `sz` per stream, and transport them to the
standard measure type `S`.

Returns a tuple `(Z, X_rest)` of the standard variates as a batch
`(getdof(μ) * prod(sz), batch dims...)` and the unconsumed rest of the
streams. The default implementation consumes variates of the size given by
[`MeasureBase.mspace_flatsize`](@ref) or
[`MeasureBase.some_mspace_elsize`](@ref), a single stream with `sz == ()`
goes through [`MeasureBase.transport_to_std_with_rest`](@ref). Measures
whose variates are composed of the variates of other measures implement
`batched_transport_to_std_with_rest` in terms of their components.
"""
function batched_transport_to_std_with_rest end

function batched_transport_to_std_with_rest(::Type{S}, μ, X::AbstractArray, sz::Dims) where {S<:StdMeasure}
    _to_std_with_rest_default(S, μ, X, sz)
end

function _to_std_with_rest_default(::Type{S}, μ, x::AbstractVector, ::Tuple{}) where {S}
    z, _, x_rest = transport_to_std_with_rest(S, μ, x)
    return z, x_rest
end
function _to_std_with_rest_default(::Type{S}, μ, X::AbstractArray, sz::Dims) where {S}
    vsz = _stream_consume_size(μ)
    X_μ, X_rest = _batched_consume(X, vsz, sz)
    Z = batched_transport_to_std(S, μ, _consumed_variates(X_μ, vsz))
    return _merge_multiplicity(Z, sz), X_rest
end
@inline _consumed_variates(X_μ::AbstractArray, ::Tuple{}) = _drop_stdstream_dim(X_μ)
@inline _consumed_variates(X_μ::AbstractArray, ::SizeLike) = X_μ

# Standard variates of `prod(sz)` variates per stream, `(dof, sz..., batch
# dims...)`, as one stream chunk `(dof * prod(sz), batch dims...)`, and
# back:
@inline _merge_multiplicity(Z::AbstractArray, sz::Dims) = _merge_leading_dims(Z, static(1) + static(length(sz)))
@inline _split_multiplicity(Z::AbstractArray, ::Tuple{}, n) = Z
@inline function _split_multiplicity(Z::AbstractArray, sz::Dims, n)
    _reshape_batch(Z, (n, sz..., Base.tail(_batch_dims(Z))...))
end


"""
    MeasureBase.batched_transport_from_std_with_rest(::Type{S}, μ, Z::AbstractArray, sz::Dims)

Consume standard variates of type `S` for a batch of variates of size
`sz` per stream from the batch `Z` of streams of standard variates (first
dimension along the streams) and transport them to `μ`.

Returns a tuple `(X, Z_rest)` of the flat batch `(flat variate dims...,
sz..., batch dims...)` of variates of `μ` and the unconsumed rest of the
streams. The default implementation consumes [`MeasureBase.fast_dof(μ)`](@ref)
entries per variate, a single stream with `sz == ()` goes through
[`MeasureBase.transport_from_std_with_rest`](@ref). Measures whose
variates are composed of the variates of other measures implement
`batched_transport_from_std_with_rest` in terms of their components.
"""
function batched_transport_from_std_with_rest end

function batched_transport_from_std_with_rest(::Type{S}, μ, Z::AbstractArray, sz::Dims) where {S<:StdMeasure}
    _from_std_with_rest_default(S, μ, Z, sz)
end

_from_std_with_rest_default(::Type{S}, μ, z::AbstractVector, ::Tuple{}) where {S} = transport_from_std_with_rest(S, μ, z)
function _from_std_with_rest_default(::Type{S}, μ, Z::AbstractArray, sz::Dims) where {S}
    _batched_from_std_bydof(S, μ, Z, sz, fast_dof(μ))
end

function _batched_from_std_bydof(::Type{S}, μ, Z::AbstractArray, sz::Dims, n::IntegerLike) where {S}
    Z_μ, Z_rest = _batched_split(Z, _chunk_rows(n, sz))
    return batched_transport_from_std(S, μ, _split_multiplicity(Z_μ, sz, n)), Z_rest
end
@noinline function _batched_from_std_bydof(::Type{S}, μ, ::AbstractArray, ::Dims, ::AbstractNoDOF) where {S}
    throw(ArgumentError("Batched transport from standard measures requires measures of type $(nameof(typeof(μ))) to have fast degrees of freedom or to implement MeasureBase.batched_transport_from_std_with_rest"))
end


"""
    MeasureBase.batched_transport_def(ν, μ, X)

Transport the flat batch `X` of variates of `μ` to a flat batch of
variates of `ν`, via the standard measure type the preferences of `ν` and
`μ` promote to. Specialize for pairs of measure types with a direct
batched transport.
"""
function batched_transport_def end

function batched_transport_def(ν, μ, X)
    S = _transport_pivot(ν, μ)
    Z = batched_transport_to_std(S, μ, X)
    Y, Z_rest = batched_transport_from_std_with_rest(S, ν, Z, ())
    if size(Z_rest, 1) != 0
        throw(ArgumentError("Degrees of freedom of source and target measure of a transport don't match"))
    end
    return Y
end


# Broadcasting a transport function over an array of variates with flat
# storage, or over the flat storage of a batch, transports the batch as a
# whole. Fused broadcast arguments are materialized first, static arrays
# are transported point by point.
function Broadcast.broadcasted(f::TransportFunction, X::AbstractArray)
    _broadcast_transport(f, X, _flat_storage(X), _static_ndims(f.μ), _static_ndims(f.ν))
end

function Broadcast.broadcasted(f::TransportFunction, bc::Broadcast.Broadcasted)
    Broadcast.broadcasted(f, Broadcast.materialize(bc))
end

# Static arrays of scalar variates are transported point by point:
Broadcast.broadcasted(f::TransportFunction, X::StaticArray) = _broadcast_static(f, X, _static_ndims(f.μ))
_broadcast_static(f::TransportFunction, X::StaticArray, ::StaticInteger{0}) = map(_Pointwise(f), X)
function _broadcast_static(f::TransportFunction, X::StaticArray, k)
    _broadcast_transport(f, X, X, k, _static_ndims(f.ν))
end

function _broadcast_transport(f::TransportFunction, X, X_flat::AbstractArray, ::StaticInteger, ::StaticInteger{K}) where {K}
    Y_flat = batched_transport_def(f.ν, f.μ, X_flat)
    return _batch_variates(Y_flat, f.ν, Val(K))
end

# Batches of tuple and named tuple variates are tuples of batches, the
# target layout follows from the target measure:
function _broadcast_transport(f::TransportFunction, X, X_flat::Union{Tuple,NamedTuple}, ::Any, ::Any)
    _structured_variates(batched_transport_def(f.ν, f.μ, X_flat), f.ν)
end
function _broadcast_transport(f::TransportFunction{<:ProductMeasure{<:Union{Tuple,NamedTuple}}}, X, X_flat::AbstractArray, ::StaticInteger, ::NoMSpaceElementSize)
    _structured_variates(batched_transport_def(f.ν, f.μ, X_flat), f.ν)
end
@inline _structured_variates(Y::Union{Tuple,NamedTuple}, ν) = _pwr_variate(ν, Y)
@inline _structured_variates(Y::AbstractArray, ν) = _batch_variates(Y, ν, Val(dynamic(_static_ndims(ν))))

_broadcast_transport(f::TransportFunction, X, ::Any, ::Any, ::Any) = map(_Pointwise(f), X)

# Prevents re-entering the broadcast hook from `map` implementations that
# broadcast (e.g. GPU arrays):
struct _Pointwise{F} <: Function
    f::F
end
@inline (p::_Pointwise)(x) = p.f(x)

# The batch of variates in the layout of the target measure over the flat
# result, nested powers included, batches of tuple variates as struct
# arrays:
@inline _batch_variates(Y::AbstractArray, ν, ::Val{K}) where {K} = _nest_batch(Y, Val(K))
@inline function _batch_variates(Y::AbstractArray, ν::PowerMeasure, ::Val)
    sliced(_pwr_variate(ν, Y), Val(length(pwr_axes(ν))))
end
@inline _nest_batch(Y::AbstractArray, ::Val{0}) = Y
@inline _nest_batch(Y::AbstractArray{<:Any,K}, ::Val{K}) where {K} = Y
@inline _nest_batch(Y::AbstractArray, ::Val{K}) where {K} = sliced(Y, Val(K))
