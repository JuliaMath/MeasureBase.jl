# A lots of this is about bridging Static and StaticArrays, both have their
# own SUnitRange and SOneTo. Also provides tools to control static vs dynamic
# array, size and axes handling.

"""
    MeasureBase.StaticUnitRange

The MeasureBase default type for static unit ranges.
"""
const StaticUnitRange = @static if isdefined(StaticArrays, :SUnitRange)
    # Unclear if StaticArrays.SUnitRange is part of StaticArrays stable API.
    # Some packages use it, but let's be careful in case it disappears.
    StaticArrays.SUnitRange
else
    Static.SUnitRange
end

"""
    MeasureBase.StaticOneTo

The MeasureBase default type for static one-based unit ranges.
"""
const StaticOneTo{T} = StaticArrays.SOneTo{T}

"""
    MeasureBase.IntegerLike

Equivalent to `Union{Integer,Static.StaticInteger}`.
"""
const IntegerLike = Union{Integer,Static.StaticInteger}

"""
    MeasureBase.SizeLike

Something that can represent the size of a collection.
"""
const SizeLike = Union{Tuple{},Tuple{Vararg{IntegerLike}},StaticArrays.Size}

"""
    MeasureBase.StaticSizeLike

Something that can represent the size of a statically sized collection.
"""
const StaticSizeLike = Union{Tuple{Vararg{StaticInteger}},StaticArrays.Size}

"""
    MeasureBase.AxesLike

Something that can represent axes of a collection.
"""
const AxesLike = Union{Tuple{},Tuple{Vararg{AbstractVector{<:IntegerLike}}}}

"""
    MeasureBase.StaticAxesLike

Something that can represent axes of a statically sized collection.
"""
@static if isdefined(StaticArrays, :SUnitRange)
    const StaticAxesLike = Union{
        Tuple{Vararg{Union{StaticArrays.SOneTo,StaticArrays.SUnitRange,Static.SUnitRange}}},
    }
else
    const StaticAxesLike =
        Union{Tuple{Vararg{Union{StaticArrays.SOneTo,Static.SUnitRange}}}}
end

"""
    const OneToLike

Alias for unit ranges that start at one.
"""
const OneToLike = Union{Base.OneTo,StaticArrays.SOneTo,Static.SOneTo}

"""
    const StaticOneToLike{N}

A static unit range from one to N.
"""
const StaticOneToLike{N} = Union{StaticArrays.SOneTo{N},Static.SOneTo{N}}

"""
    const StaticUnitRangeLike

A static unit range.
"""
@static if isdefined(StaticArrays, :SUnitRange)
    const StaticUnitRangeLike =
        Union{StaticArrays.SOneTo,StaticArrays.SUnitRange,Static.SUnitRange}
else
    const StaticUnitRangeLike = Union{StaticArrays.SOneTo,Static.SUnitRange}
end

"""
    MeasureBase.one_to(n::IntegerLike)

Creates a range from one to n.

Returns an instance of `Base.OneTo` or `Static.SOneTo`, depending
on the type of `n`.
"""
@inline one_to(n::Integer) = Base.OneTo(n)
@inline one_to(::Static.StaticInteger{N}) where {N} = Static.SOneTo{N}()

"""
    MeasureBase.asnonstatic(x)

Return a non-static equivalent of `x`.

Defaults to `Static.dynamic(x)`.
"""
@inline asnonstatic(x::Number) = dynamic(x)
@inline asnonstatic(::Tuple{}) = ()
@static if isdefined(StaticArrays, :SUnitRange)
    @inline asnonstatic(r::StaticArrays.SUnitRange) = r[begin]:r[end]
end
@inline asnonstatic(r::AbstractUnitRange) = asnonstatic(r[begin]):asnonstatic(r[end])
@inline asnonstatic(r::Base.OneTo) = Base.OneTo(asnonstatic(r.stop))
@inline asnonstatic(::StaticOneToLike{N}) where {N} = Base.OneTo(N)
@inline asnonstatic(x::SizeLike) = map(asnonstatic, x)
@inline asnonstatic(::StaticArrays.Size{TPL}) where {TPL} = TPL
@inline asnonstatic(x::AxesLike) = map(asnonstatic, x)

"""
    MeasureBase.fill_with(x, sz::NTuple{N,<:IntegerLike}) where N

Creates an array of size `sz` filled with `x`.

The result will typically be either a `FillArrays.Fill` or a static array,
"""
function fill_with end

@inline fill_with(x::T, n::IntegerLike) where {T} = fill_with(x, (n,))

@inline fill_with(x::T, ::Tuple{}) where {T} = FillArrays.Fill(x)

@inline fill_with(x, sz::SizeLike) = fill_with(x, size2axes(sz))

@inline function fill_with(x::T, sz::StaticSizeLike) where {T}
    fill(x, staticarray_type(T, canonical_size(sz)))
end

@inline function fill_with(x, axs::AxesLike)
    dyn_axs = map(asnonstatic, axs)
    FillArrays.Fill(x, dyn_axs)
end

# While `FillArrays.Fill` (mostly?) works with axes that are static unit
# ranges, some operations that automatic differentiation requires do fail
# on such instances of `Fill` (e.g. `reshape` from dynamic to static size).
# So need to build a filled static array:
@inline function fill_with(x::T, axs::Tuple{Vararg{StaticOneToLike}}) where {T}
    sz = axes2size(axs)
    fill(x, staticarray_type(T, sz))
end

"""
    MeasureBase.staticarray_type(T, sz::StaticArrays.Size)

Returns the type of a static array with element type `T` and size `sz`.
"""
function staticarray_type end

@inline @generated function staticarray_type(
    ::Type{T},
    ::StaticArrays.Size{sz},
) where {T,sz}
    N = length(sz)
    len = prod(sz)
    :(SArray{Tuple{$sz...},T,$N,$len})
end

"""
    MeasureBase.maybestatic_reshape(A, sz)

Reshapes array `A` to sizes `sz`.

If `A` is a static array and `sz` is static, the result is a static array.
"""
function maybestatic_reshape end

maybestatic_reshape(A, sz) = reshape(A, canonical_size(sz))
function maybestatic_reshape(A, sz::StaticSizeLike)
    StaticArrays.SArray(reshape(A, canonical_size(sz)))
end
function maybestatic_reshape(A::StaticArray, sz::Tuple{Vararg{StaticInteger}})
    staticarray_type(eltype(A), canonical_size(sz))(Tuple(A))
end

"""
    MeasureBase.maybestatic_length(x)

Returns the length of `x` as a dynamic or static integer.
"""
@inline maybestatic_length(::Number) = static(1)
@inline maybestatic_length(::Tuple{}) = static(0)
@inline maybestatic_length(::Tuple{Vararg{Any,N}}) where {N} = static(N)
@inline maybestatic_length(nt::NamedTuple) = maybestatic_length(values(nt))
@inline maybestatic_length(A::AbstractArray) = size2length(maybestatic_size(A))
@static if isdefined(StaticArrays, :SUnitRange)
    @inline maybestatic_length(r::StaticArrays.SUnitRange) =
        maybestatic_last(r) - maybestatic_first(r) + static(1)
end
@inline maybestatic_length(r::AbstractUnitRange) =
    maybestatic_last(r) - maybestatic_first(r) + static(1)
@inline maybestatic_length(r::Base.OneTo) = length(r)
@inline maybestatic_length(::StaticArrays.SOneTo{N}) where {N} = static(N)
@inline maybestatic_length(::Static.SOneTo{N}) where {N} = static(N)

"""

    MeasureBase.maybestatic_size(x)

Returns the size of `x` as a tuple of dynamic or static integers.
"""
@inline maybestatic_size(::Number) = ()
@inline maybestatic_size(::Tuple{}) = StaticArrays.Size(0)
@inline maybestatic_size(::Tuple{Vararg{Any,N}}) where {N} = StaticArrays.Size(N)
@inline maybestatic_size(nt::NamedTuple) = maybestatic_size(values(nt))
@inline maybestatic_size(::StaticArrays.Size{tpl}) where {tpl} =
    StaticArrays.Size(length(tpl))
@inline maybestatic_size(A::AbstractArray) = axes2size(maybestatic_axes(A))
@inline maybestatic_size(A::StaticArray) = StaticArrays.Size(A)

"""
    MeasureBase.maybestatic_axes(x)::Tuple{Vararg{IntegerLike}}

Returns the size of `x` as a tuple of dynamic or static integers.
"""
@inline maybestatic_axes(::Number) = ()

@inline maybestatic_axes(::Tuple{}) = (StaticOneTo(0),)
@inline maybestatic_axes(::Tuple{Vararg{Any,N}}) where {N} = (StaticOneTo(N),)
@inline maybestatic_axes(nt::NamedTuple) = maybestatic_axes(values(nt))
@inline maybestatic_axes(::StaticArrays.Size{tpl}) where {tpl} = (StaticOneTo(length(tpl)),)
@inline maybestatic_axes(::StaticOneToLike{N}) where {N} = (StaticOneTo(N),)
@static if isdefined(StaticArrays, :SUnitRange)
    @inline maybestatic_axes(r::StaticArrays.SUnitRange) = axes(r)
end
@inline maybestatic_axes(r::Static.OptionallyStaticUnitRange) = canonical_axes(axes(r))
@inline maybestatic_axes(r::AbstractUnitRange) = axes(r)
@inline maybestatic_axes(A::AbstractArray) = axes(A)
@inline maybestatic_axes(A::StaticArray) = axes(A)

"""
    MeasureBase.axes2size(x::Tuple)
    MeasureBase.axes2size(x::StaticArrays.Size)

Get a length from a size (tuple).
"""
@inline axes2size(::Tuple{}) = ()
@inline axes2size(axs::Tuple) = canonical_size(map(maybestatic_length, axs))

"""map(maybestatic_length, axs)
    MeasureBase.size2axes(sz::Tuple)
    MeasureBase.size2axes(sz::StaticArrays.Size)

Get one-based indexing axes from a size.
"""
@inline size2axes(::Tuple{}) = ()
@inline size2axes(sz::Tuple) = canonical_axes(map(one_to, sz))
@inline size2axes(::StaticArrays.Size{TPL}) where {TPL} = map(StaticOneTo, TPL)

"""
    MeasureBase.size2length(sz::Tuple)
    MeasureBase.size2length(sz::StaticArrays.Size)

Get a length from a size (tuple).
"""
@inline size2length(::Tuple{}) = static(1)
@inline size2length(sz::Tuple) = prod(sz)
@inline size2length(::StaticArrays.Size{TPL}) where {TPL} = static(prod(TPL))

"""
    MeasureBase.asaxes(axs::AxesLike)
    MeasureBase.asaxes(sz::SizeLike)
    MeasureBase.asaxes(len::IntegerLike)

Converts axes or a size or a length of a collection to axes.

One-based indexing will be used if the indexing offset can't be inferred from
the given dimensions.
"""
@inline asaxes(::Tuple{}) = ()
@inline asaxes(axs::AxesLike) = axs
@inline asaxes(sz::SizeLike) = size2axes(sz)
@inline asaxes(len::IntegerLike) = size2axes((len,))

"""
    MeasureBase.maybestatic_eachindex(x)

Returns the the index range of `x` as a dynamic or static integer range
"""
maybestatic_eachindex(::Tuple{}) = StaticOneTo(0)
maybestatic_eachindex(::Tuple{Vararg{Any,N}}) where {N} = StaticOneTo(N)
maybestatic_eachindex(nt::NamedTuple) = maybestatic_eachindex(values(nt))
maybestatic_eachindex(x::AbstractArray) = canonical_indices(eachindex(x))

"""
    MeasureBase.maybestatic_first(A)

Returns the first element of `A` as a dynamic or static value.
"""
maybestatic_first(tpl::Tuple) = tpl[begin]
maybestatic_first(nt::NamedTuple) = nt[begin]
maybestatic_first(A::AbstractArray) = A[begin]
maybestatic_first(::StaticArrays.Size{tpl}) where {tpl} = static(tpl[begin])
maybestatic_first(::StaticArrays.SOneTo{N}) where {N} = static(1)
@static if isdefined(StaticArrays, :SUnitRange)
    maybestatic_first(::StaticArrays.SUnitRange{B,L}) where {B,L} = static(B)
end
function maybestatic_first(
    ::Static.OptionallyStaticUnitRange{<:Static.StaticInteger{from},<:Static.StaticInteger},
) where {from}
    static(from)
end

"""
    MeasureBase.maybestatic_last(A)

Returns the last element of `A` as a dynamic or static value.
"""
maybestatic_last(tpl::Tuple) = tpl[end]
maybestatic_last(nt::NamedTuple) = nt[end]
maybestatic_last(A::AbstractArray) = A[end]
maybestatic_last(::StaticArrays.Size{tpl}) where {tpl} = static(tpl[end])
maybestatic_last(::StaticArrays.SOneTo{N}) where {N} = static(N)
@static if isdefined(StaticArrays, :SUnitRange)
    maybestatic_last(::StaticArrays.SUnitRange{B,L}) where {B,L} = static(B + L - 1)
end
function maybestatic_last(
    ::Static.OptionallyStaticUnitRange{<:Any,<:Static.StaticInteger{until}},
) where {until}
    static(until)
end

"""
    MeasureBase.canonical_indices(idxs::AbstractVector{<:IntegerLike})

Return the canonical representation of a collection axis indices.
"""
@inline canonical_indices(idxs::AbstractVector{<:IntegerLike}) = idxs
@inline canonical_indices(idxs::AbstractArray{<:CartesianIndex}) = idxs
@inline canonical_indices(
    ::Static.OptionallyStaticUnitRange{<:StaticInteger{1},<:StaticInteger{N}},
) where {N} = StaticArrays.SOneTo{N}()
@inline canonical_indices(
    ::Static.OptionallyStaticUnitRange{<:StaticInteger{A},<:StaticInteger{B}},
) where {A,B} = StaticUnitRange(A, B)
@inline canonical_indices(
    r::Static.OptionallyStaticUnitRange{<:StaticInteger{1},<:Integer},
) = Base.OneTo(last(r))

"""
    MeasureBase.canonical_size(sz::SizeLike)

Return the canonical representation of a collection size.
"""
@inline canonical_size(sz::SizeLike) = sz
@inline canonical_size(sz::Tuple{Vararg{Static.StaticInteger}}) =
    StaticArrays.Size{map(dynamic, sz)}()

"""
    MeasureBase.canonical_axes(sz::SizeLike)

Return the canonical representation collection axes.
"""
@inline canonical_axes(axs::AxesLike) = map(canonical_indices, axs)
