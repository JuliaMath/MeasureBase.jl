"""
    MeasureBase.IntegerLike

Equivalent to `Union{Integer,Static.StaticInt}`.
"""
const IntegerLike = Union{Integer,Static.StaticInt}


"""
    MeasureBase.one_to(n::IntegerLike)

Creates a range from one to n.

Returns an instance of `Base.OneTo` or `Static.SOneTo`, depending
on the type of `n`.
"""
@inline one_to(n::Integer) = Base.OneTo(n)
@inline one_to(::Static.StaticInt{N}) where N = Static.SOneTo{N}()


_dynamic(x::Number) = dynamic(x)
_dynamic(::Static.SOneTo{N}) where N = Base.OneTo(N)
_dynamic(r::AbstractUnitRange) = minimum(r):maximum(r)


"""
    MeasureBase.fill_with(x, sz::NTuple{N,<:IntegerLike}) where N

Creates an array of size `sz` filled with `x`.

Returns an instance of `FillArrays.Fill`.
"""
function fill_with end

@inline fill_with(x::T, sz::Tuple{Vararg{<:IntegerLike,N}}) where {T,N} = fill_with(x, map(one_to, sz))

@inline function fill_with(x::T, axs::Tuple{Vararg{<:AbstractUnitRange,N}}) where {T,N}
    # While `FillArrays.Fill` (mostly?) works with axes that are static unit
    # ranges, some operations that automatic differentiation requires do fail
    # on such instances of `Fill` (e.g. `reshape` from dynamic to static size).
    # So need to use standard ranges for the axes for now:
    dyn_axs = map(_dynamic, axs)
    FillArrays.Fill(x, dyn_axs)
end


"""
    MeasureBase.dslength(x)::IntegerLike

Returns the length of `x` as a dynamic or static integer.
"""
dslength(x) = length(x)
dslength(x::AbstractUnitRange) = length(x)
dslength(::Static.OptionallyStaticUnitRange{StaticInt{A}, StaticInt{B}}) where {A,B} = StaticInt{B - A + 1}()


"""
    MeasureBase.dssize(x)::Tuple{Vararg{IntegerLike}}

Returns the size of `x` as a tuple of dynamic or static integers.
"""
dssize(x) = size(x)
