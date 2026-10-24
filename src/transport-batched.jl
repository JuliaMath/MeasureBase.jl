# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Batched-first transport over flat batches of variates `(variate dims...,
# batch dims...)`, zero batch dims meaning a single variate. Streams of
# standard variates are batches `(dof, batch dims...)`, consumed along their
# first dimension. Kernels know the stream rank of their measure (see
# `stream_ndims`), structural measures implement them once in terms of the
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
    _batched_to_std(S, μ, X, _stream_ndims(μ))
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

# A flat batch of variates of `μ` as a batch of streams, the variate
# dimensions merged into the first dimension. Tuples of batches (tuple
# products and their powers) interleave the rows of their components
# variate by variate.
@inline _as_stream_batch(X, μ) = _as_stream_batch(X, _stream_ndims(μ))
@inline function _as_stream_batch(X::AbstractArray, ::StaticInteger{K}) where {K}
    ndims(X) >= K || _throw_size_mismatch()
    merge_leading_dims(X, static(K))
end
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
    n = prod(asnonstatic(_pwr_dims(μ)))
    parts = map(values(X), values(marginals(ν))) do Xi, m
        A = _as_stream_batch(Xi, m)
        reshape(A, (size(A, 1), n, ntuple(i -> size(A, 1 + n_pwr + i), Val(ndims(A) - 1 - n_pwr))...))
    end
    merge_leading_dims(vcat(parts...), static(2))
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
    _batched_from_std(S, μ, Z, _stream_ndims(μ))
end

@inline _batched_from_std(::Type{S}, μ, Z::AbstractArray, ::StaticInteger{0}) where {S} = _from_std_scalar(S, μ, Z)
@inline _batched_from_std(::Type{S}, μ, Z::AbstractArray, ::StaticInteger{K}) where {S,K} = _from_std_columns(S, μ, Z)
@noinline function _batched_from_std(::Type{S}, μ, ::AbstractArray, ::NoMSpaceElementSize) where {S}
    throw(ArgumentError("Batched transport requires MeasureBase.mspace_ndims to be declared for measures of type $(nameof(typeof(μ))) or MeasureBase.batched_transport_from_std to be implemented"))
end

@inline _from_std_scalar(::Type{S}, μ, z::AbstractVector) where {S} = transport_from_std(S, μ, _single_std_element(z))
@inline function _from_std_scalar(::Type{S}, μ, Z::AbstractArray) where {S}
    broadcast(Base.Fix1(_FromStd{S}(), μ), _drop_stdstream_dim(Z))
end
@inline _from_std_columns(::Type{S}, μ, z::AbstractVector) where {S} = transport_from_std(S, μ, z)
@inline function _from_std_columns(::Type{S}, μ, Z::AbstractArray) where {S}
    stacked(map(Base.Fix1(_FromStd{S}(), μ), sliced(Z, Val(1))))
end


"""
    MeasureBase.batched_transport_to_std_with_rest(::Type{S}, μ, X::AbstractArray, sz::SizeLike)

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

function batched_transport_to_std_with_rest(::Type{S}, μ, X::AbstractArray, sz::SizeLike) where {S<:StdMeasure}
    _to_std_with_rest_default(S, μ, X, sz)
end

function _to_std_with_rest_default(::Type{S}, μ, x::AbstractVector, ::Tuple{}) where {S}
    z, _, x_rest = transport_to_std_with_rest(S, μ, x)
    return z, x_rest
end
function _to_std_with_rest_default(::Type{S}, μ, X::AbstractArray, sz::SizeLike) where {S}
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
@inline _merge_multiplicity(Z::AbstractArray, sz::SizeLike) = merge_leading_dims(Z, static(1) + maybestatic_length(size_tuple(sz)))
@inline _split_multiplicity(Z::AbstractArray, ::Tuple{}, n) = Z
@inline function _split_multiplicity(Z::AbstractArray, sz::SizeLike, n)
    maybestatic_reshape(Z, (n, size_tuple(sz)..., Base.tail(_batch_dims(Z))...))
end


"""
    MeasureBase.batched_transport_from_std_with_rest(::Type{S}, μ, Z::AbstractArray, sz::SizeLike)

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

function batched_transport_from_std_with_rest(::Type{S}, μ, Z::AbstractArray, sz::SizeLike) where {S<:StdMeasure}
    _from_std_with_rest_default(S, μ, Z, sz)
end

_from_std_with_rest_default(::Type{S}, μ, z::AbstractVector, ::Tuple{}) where {S} = transport_from_std_with_rest(S, μ, z)
function _from_std_with_rest_default(::Type{S}, μ, Z::AbstractArray, sz::SizeLike) where {S}
    _batched_from_std_bydof(S, μ, Z, sz, fast_dof(μ))
end

function _batched_from_std_bydof(::Type{S}, μ, Z::AbstractArray, sz::SizeLike, n::IntegerLike) where {S}
    Z_μ, Z_rest = _batched_split(Z, _chunk_rows(n, sz))
    return batched_transport_from_std(S, μ, _split_multiplicity(Z_μ, sz, n)), Z_rest
end
# Without fast degrees of freedom, the exact ones split the streams:
function _batched_from_std_bydof(::Type{S}, μ, Z::AbstractArray, sz::SizeLike, ::NoFastDOF) where {S}
    _batched_from_std_bydof(S, μ, Z, sz, getdof(μ))
end
function _from_std_with_rest_bydof(::Type{S}, μ, z::AbstractVector, ::NoFastDOF) where {S}
    _from_std_with_rest_bydof(S, μ, z, getdof(μ))
end
@noinline function _batched_from_std_bydof(::Type{S}, μ, ::AbstractArray, ::SizeLike, ::NoDOF) where {S}
    throw(ArgumentError("Batched transport from standard measures requires measures of type $(nameof(typeof(μ))) to have degrees of freedom or to implement MeasureBase.batched_transport_from_std_with_rest"))
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


"""
    transports(f::TransportFunction, X)

Transport each variate in the batch `X` with `f = transport_to(ν, μ)`.

`X` is a batch of variates of `μ` as for [`logdensities`](@ref): an array
of variates, the flat storage of a batch (an array of numbers with the
variate dimensions leading, see [`MeasureBase.stream_ndims`](@ref)), or a
tuple resp. `NamedTuple` of batches for tuple resp. `NamedTuple`
variates. Returns the array of the transported variates over the batch
dimensions of `X` (a single variate for a single variate `X`),
semantically equivalent to `f.(xs)` for the array `xs` of the variates in
`X`. Array variates are returned as nested arrays over flat storage
(ArraysOfArrays), tuple variates as struct arrays. Variates of
value-dependent size (e.g. of binds) are transported one by one, to an
array of variates.

Batches with flat storage are transported in fused operations,
compatible with GPU and traced arrays. The result is materialized, it is
no view of temporary arrays and shares no memory with `X` (also for
identity transports), static batches give static results.
Broadcasting `f`, `f.(X)`, is equivalent to `transports(f, X)`.
"""
function transports end
export transports

@inline transports(f::TransportFunction, X) = _transports(f, X, _flat_storage(X), _stream_ndims(f.μ), _stream_ndims(f.ν))
@inline transports(f::TransportFunction, x::Number) = f(x)
@inline transports(f::TransportFunction, X::Union{Tuple,NamedTuple}) = _transported_batch(f, X, X)
# Static arrays of scalar variates are transported point by point:
@inline transports(f::TransportFunction, X::StaticArray) = _transports_static(f, X, _stream_ndims(f.μ))
@inline _transports_static(f::TransportFunction, X::StaticArray, ::StaticInteger{0}) = map(_Pointwise(f), X)
@inline function _transports_static(f::TransportFunction, X::StaticArray, k)
    _transports(f, X, _flat_storage(X), k, _stream_ndims(f.ν))
end

function _transports(f::TransportFunction, X, X_flat::AbstractArray, ::StaticInteger, ::StaticInteger)
    _check_flatsize(X_flat, mspace_flatsize(f.μ))
    _transported_batch(f, X, X_flat)
end
# Batches of tuple and named tuple variates are tuples of batches:
function _transports(f::TransportFunction, X, X_flat::Union{Tuple,NamedTuple}, ::Any, ::Any)
    _transported_batch(f, X, X_flat)
end
_transports(f::TransportFunction, X, ::Any, k, ::Any) = _map_variates(f, X, k)

_transported_batch(f::TransportFunction, X, X_flat) = _transported_batch(f, X, X_flat, _flat_batchable(f.ν))
function _transported_batch(f::TransportFunction, X, X_flat, ::True)
    Y = _batch_variates(_materialized(batched_transport_def(f.ν, f.μ, X_flat), X_flat), f.ν)
    _check_transported_shape(Y, X)
end
# Variates of value-dependent size (of binds and measures composed of them)
# have no flat batches, they are transported from the standard variates of
# the batch one by one, to an array of variates as `rand` generates them for
# powers of such measures:
function _transported_batch(f::TransportFunction, X, X_flat, ::False)
    S = _transport_pivot(f.ν, f.μ)
    _from_std_variates(S, f.ν, batched_transport_to_std(S, f.μ, X_flat))
end
@inline _from_std_variates(::Type{S}, ν, z::AbstractVector) where {S} = _from_std_stream(S, ν, z, _static_ndims(ν))
function _from_std_variates(::Type{S}, ν, Z::AbstractArray) where {S}
    map(z -> _from_std_stream(S, ν, z, _static_ndims(ν)), sliced(Z, Val(1)))
end

# Measures without fixed stream sizes have variates of value-dependent size
# if they have no degrees of freedom either:
@inline _flat_batchable(ν) = _flat_batchable(ν, fixed_stream_size(ν))
@inline _flat_batchable(::Any, ::True) = static(true)
@inline _flat_batchable(ν, ::False) = _has_dof(getdof(ν))
@inline _has_dof(::Any) = static(true)
@inline _has_dof(::NoDOF) = static(false)

# Arrays of variates must keep their shape, a mismatch means the variate
# rank of the measure doesn't match the batch:
@inline _check_transported_shape(Y, ::Union{AbstractArray{<:Number},Tuple,NamedTuple}) = Y
@inline function _check_transported_shape(Y, X::AbstractArray)
    size(Y) == size(X) || _throw_size_mismatch()
    return Y
end

Broadcast.broadcasted(f::TransportFunction, X::AbstractArray) = transports(f, X)

# Fused broadcast arguments are materialized once, arrays of fixed-size
# array variates directly into flat storage:
function Broadcast.broadcasted(f::TransportFunction, bc::Broadcast.Broadcasted)
    transports(f, _materialized_arg(bc, fixed_stream_size(f.μ)))
end
@inline _materialized_arg(bc::Broadcast.Broadcasted, ::True) = _materialized_arg(bc, Broadcast.combine_eltypes(bc.f, bc.args))
@inline _materialized_arg(bc::Broadcast.Broadcasted{<:StaticArrays.StaticArrayStyle}, ::True) = Broadcast.materialize(bc)
@inline _materialized_arg(bc::Broadcast.Broadcasted, ::Any) = Broadcast.materialize(bc)
@inline _materialized_arg(bc::Broadcast.Broadcasted, ::Type{<:AbstractArray{<:Number}}) = stack(bc)

# Flat batches of array variates are mapped over their variate slices:
_map_variates(f::TransportFunction, X, ::Any) = map(_Pointwise(f), X)
function _map_variates(f::TransportFunction, X::AbstractArray{<:Number,N}, ::StaticInteger{K}) where {N,K}
    N >= K || _throw_size_mismatch()
    if K == 0
        map(_Pointwise(f), X)
    elseif N == K
        f(X)
    else
        map(_Pointwise(f), sliced(X, Val(K)))
    end
end

# Prevents re-entering the broadcast hook from `map` implementations that
# broadcast (e.g. GPU arrays):
struct _Pointwise{F} <: Function
    f::F
end
@inline (p::_Pointwise)(x) = p.f(x)

# The batch of variates of `ν` in their layout over the flat batch `Y`,
# nested powers included, batches of tuple variates as struct arrays.
# Tuples of batches of other measures have the layout of their base
# measure, which shares their variate space:
@inline _batch_variates(Y, ν) = _batch_variates(Y, ν, _static_ndims(ν))
@inline _batch_variates(Y, ν::PowerMeasure) = _nest_batch(_pwr_variate(ν, Y), Val(length(pwr_axes(ν))))
@inline _batch_variates(Y::Union{Tuple,NamedTuple}, ν::ProductMeasure{<:Union{Tuple,NamedTuple}}) = _single_or_batch(_pwr_variate(ν, Y))
@inline _batch_variates(Y::AbstractArray, ν, ::StaticInteger{K}) where {K} = _nest_batch(Y, Val(K))
@inline _batch_variates(Y::Union{Tuple,NamedTuple}, ν, ::Any) = _batch_variates_via(basemeasure(ν), ν, Y)
@inline _batch_variates_via(β, ::Any, Y) = _batch_variates(Y, β)
@noinline function _batch_variates_via(::M, ::M, Y) where {M}
    throw(ArgumentError("Batches of tuple variates of measures of type $(nameof(M)) are not supported"))
end

# A single variate (no batch dimensions) is returned as such:
@inline _single_or_batch(A::AbstractArray{<:Any,0}) = A[]
@inline _single_or_batch(A::AbstractArray) = A

@inline _nest_batch(Y::AbstractArray, ::Val{0}) = Y
@inline _nest_batch(Y::AbstractArray{<:Any,K}, ::Val{K}) where {K} = Y
@inline _nest_batch(Y::AbstractArray, ::Val{K}) where {K} = sliced(Y, Val(K))


"""
    transports_with_ladj(f::TransportFunction, X)

Transport each variate in the batch `X` with `f = transport_to(ν, μ)`
and compute the log-abs-det Jacobians of the transport, the batched form
of `ChangesOfVariables.with_logabsdet_jacobian(f, x)`.

Returns a tuple `(Y, ladjs)` of the transported batch `Y`, as returned by
[`transports(f, X)`](@ref transports), and the array `ladjs` of the
log-abs-det Jacobians over the batch dimensions of `X` (a number for a
single variate), like the results of [`logdensities`](@ref). The
Jacobians are the ratios of the normalized densities of `μ` and `ν` (see
[`transport_to`](@ref)), evaluated in batched density evaluations.
"""
function transports_with_ladj end
export transports_with_ladj

function transports_with_ladj(f::TransportFunction, X)
    Y = transports(f, X)
    ladjs = broadcast(_transport_ladj, logdensities(_unweighted(f.μ), X), logdensities(_unweighted(f.ν), Y))
    return Y, ladjs
end
