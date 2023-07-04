"""
    MeasureBase.IntegerLike

Equivalent to `Union{Integer,Static.StaticInteger}`.
"""
const IntegerLike = Union{Integer,Static.StaticInteger}


"""
    const UnitRangeFromOne

Alias for unit ranges that start at one.
"""
const UnitRangeFromOne = Union{Base.OneTo, Static.OptionallyStaticUnitRange, StaticArrays.SOneTo}


"""
    const StaticOneTo{N}

A static unit range from one to N.
"""
const StaticOneTo{N} = Union{Static.OptionallyStaticUnitRange{StaticInt{1},StaticInt{N}}, StaticArrays.SOneTo{N}}


"""
    const StaticUnitRange

A static unit range.
"""
const StaticUnitRange = Union{Static.OptionallyStaticUnitRange{<:StaticInt,<:StaticInt}, StaticArrays.SOneTo}


"""
    MeasureBase.one_to(n::IntegerLike)

Creates a range from one to n.

Returns an instance of `Base.OneTo` or `Static.SOneTo`, depending
on the type of `n`.
"""
@inline one_to(n::Integer) = Base.OneTo(n)
@inline one_to(::Static.StaticInteger{N}) where {N} = Static.SOneTo{N}()

_dynamic(x::Number) = dynamic(x)
_dynamic(::Static.SOneTo{N}) where {N} = Base.OneTo(N)
_dynamic(r::AbstractUnitRange) = minimum(r):maximum(r)

"""
    MeasureBase.maybestatic_fill(x, sz::NTuple{N,<:IntegerLike}) where N

Creates an array of size `sz` filled with `x`.

Returns an instance of `FillArrays.Fill`.
"""
function maybestatic_fill end

@inline maybestatic_fill(x::T, ::Tuple{}) where T = FillArrays.Fill(x)

@inline function maybestatic_fill(x::T, sz::Tuple{Vararg{IntegerLike,N}}) where {T,N}
    maybestatic_fill(x, map(one_to, sz))
end

@inline function maybestatic_fill(x::T, axs::Tuple{Vararg{AbstractUnitRange,N}}) where {T,N}
    # While `FillArrays.Fill` (mostly?) works with axes that are static unit
    # ranges, some operations that automatic differentiation requires do fail
    # on such instances of `Fill` (e.g. `reshape` from dynamic to static size).
    # So need to use standard ranges for the axes for now:
    dyn_axs = map(_dynamic, axs)
    FillArrays.Fill(x, dyn_axs)
end

@inline function  maybestatic_fill(x::T, axs::Tuple{Vararg{StaticOneTo}}) where T
    fill(x, staticarray_type(T, map(maybestatic_length, axs)))
end

@inline function  maybestatic_fill(x::T, sz::Tuple{Vararg{StaticInteger}}) where T
    fill(x, staticarray_type(T, sz))
end


"""
    staticarray_type(T, sz::Tuple{Vararg{StaticInteger}})

Returns the type of a static array with element type `T` and size `sz`.
"""
function staticarray_type end

@inline @generated function staticarray_type(::Type{T}, sz::Tuple{Vararg{StaticInteger,N}}) where {T,N}
    szs = map(p -> p.parameters[1], sz.parameters)
    len = prod(szs)
    :(SArray{Tuple{$szs...},T,$N,$len})
end


"""
    MeasureBase.maybestatic_reshape(A, sz)

Reshapes array `A` to sizes `sz`.

If `A` is a static array and `sz` is static, the result is a static array.
"""
function maybestatic_reshape end

maybestatic_reshape(A, sz) = reshape(A, sz)
maybestatic_reshape(A::StaticArray, sz::Tuple{Vararg{StaticInteger}}) = _sarray_type(eltype(A), sz)(Tuple(A))


"""
    MeasureBase.maybestatic_length(x)::IntegerLike

Returns the length of `x` as a dynamic or static integer.
"""
maybestatic_length(x) = length(x)
maybestatic_length(x::AbstractUnitRange) = length(x)
maybestatic_length(::Tuple{Vararg{Any,N}}) where N = static(N)
maybestatic_length(nt::NamedTuple) = maybestatic_length(values(nt))
maybestatic_length(x::StaticArray) = maybestatic_length(maybestatic_eachindex(x))
maybestatic_length(::StaticArrays.SOneTo{N}) where {N} = static(N)
function maybestatic_length(
    ::Static.OptionallyStaticUnitRange{<:StaticInteger{A},<:StaticInteger{B}},
) where {A,B}
    StaticInt{B - A + 1}()
end


"""
    MeasureBase.maybestatic_size(x)::Tuple{Vararg{IntegerLike}}

Returns the size of `x` as a tuple of dynamic or static integers.
"""
maybestatic_size(x) = map(maybestatic_length, axes(x))


"""
    MeasureBase.maybestatic_eachindex(x)

Returns the the index range of `x` as a dynamic or static integer range
"""
maybestatic_eachindex(x::AbstractArray) = _conv_static_eachindex(eachindex(x))
maybestatic_eachindex(::Tuple{Vararg{Any,N}}) where N = static(1):static(N)
maybestatic_eachindex(nt::NamedTuple) = maybestatic_eachindex(values(nt))

_conv_static_eachindex(idxs) = idxs
_conv_static_eachindex(::Static.SOneTo{N}) where {N} = static(1):static(N)


"""
    MeasureBase.maybestatic_first(A)

Returns the first element of `A` as a dynamic or static value.
"""
maybestatic_first(A::AbstractArray) = first(A)
maybestatic_first(::StaticArrays.SOneTo{N}) where N = static(1)
maybestatic_first(::Static.OptionallyStaticUnitRange{<:Static.StaticInteger{from},<:Static.StaticInteger}) where from = static(from)


"""
    MeasureBase.maybestatic_last(A)

Returns the last element of `A` as a dynamic or static value.
"""
maybestatic_last(A::AbstractArray) = last(A)
maybestatic_last(::StaticArrays.SOneTo{N}) where N = static(N)
maybestatic_last(::Static.OptionallyStaticUnitRange{<:Static.StaticInteger,<:Static.StaticInteger{until}}) where until = static(until)
