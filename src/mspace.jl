"""
    MeasureBase.NoMSpaceElementSize{MU}

Indicates that either the measurable space of measures of type `MU` is not
a space over arrays, or that the size of the arrays is not fixed or can not
be easily/efficiently determined.
"""
struct NoMSpaceElementSize{MU} end


"""
    mspace_elsize(μ)

For a measure `μ` over an array-valued measurable space, return the size of
the arrays that are the elements of the space, `()` for scalar variates.

The size is static where it is known statically. Returns
[`NoMSpaceElementSize{typeof(μ)}()`](@ref) if the elements of the space
are not arrays of one common size, e.g. for structured variates or variates
whose size depends on the value, or if the size can not be determined
efficiently.

See also [`MeasureBase.mspace_flatsize`](@ref).
"""
function mspace_elsize end
export mspace_elsize

@inline mspace_elsize(μ::AbstractMeasure) = NoMSpaceElementSize{typeof(μ)}()


"""
    MeasureBase.mspace_flatsize(μ)

Return the size of the flat storage of a variate of `μ`, `()` for scalar
variates.

Variates of powers of measures with array-valued variates are nested
arrays, their flat storage has the size of the inner arrays followed by
the size of the power. Returns [`NoMSpaceElementSize{typeof(μ)}()`](@ref)
if the variates of `μ` have no flat storage of a common size.

See also [`mspace_elsize`](@ref).
"""
function mspace_flatsize end

@inline mspace_flatsize(μ::AbstractMeasure) = NoMSpaceElementSize{typeof(μ)}()

@inline _cat_sizes(a::SizeLike, b::SizeLike) = canonical_size((_size_dims(a)..., _size_dims(b)...))
@inline _size_dims(sz::Tuple) = sz
@inline _size_dims(::StaticArrays.Size{S}) where {S} = map(static, S)
@inline _cat_sizes(a::NoMSpaceElementSize, ::SizeLike) = a
@inline _cat_sizes(::SizeLike, b::NoMSpaceElementSize) = b
@inline _cat_sizes(a::NoMSpaceElementSize, ::NoMSpaceElementSize) = a


"""
    MeasureBase.some_mspace_elsize(μ::AbstractMeasure)

For a measure `μ` over an array-valued measurable space, return the size of
an arbitrary element of the space.

Use with caution, the space of some measures is made up of arrays of
different sizes!

In general, use [`mspace_elsize(μ)`](@ref) instead. `some_mspace_elsize` is
useful if the measurable space is expected to contain only arrays of the
same size but there is no way to prove this automatically. Algorithms that
use the returned size should always check that it matches the size of each
point of the space that is processed.
"""
function some_mspace_elsize end

@inline some_mspace_elsize(μ) = _mspace_some_elsize_impl(μ, mspace_elsize(μ))

@inline _mspace_some_elsize_impl(::AbstractMeasure, sz::SizeLike) = sz
_mspace_some_elsize_impl(μ::AbstractMeasure, ::NoMSpaceElementSize) =
    maybestatic_size(testvalue(μ))

@inline _value_elsize(::Number) = ()
@inline _value_elsize(x::AbstractArray) = maybestatic_size(x)
@inline _value_elsize(x) = NoMSpaceElementSize{typeof(x)}()

@inline _value_flatsize(::Number) = ()
@inline _value_flatsize(x::AbstractArray{<:Number}) = maybestatic_size(x)
@inline _value_flatsize(x) = NoMSpaceElementSize{typeof(x)}()

@inline _scalar_or_unknown(::Tuple{}) = ()
@inline _scalar_or_unknown(sz::NoMSpaceElementSize) = sz
@inline _scalar_or_unknown(sz) = NoMSpaceElementSize{typeof(sz)}()


"""
    MeasureBase.mspace_flatsize(::Type{MU})

The flat variate size of measures of type `MU`, if it is determined by the
type alone, e.g. `()` for measures with scalar variates. Returns
`NoMSpaceElementSize{MU}()` otherwise.

Composite measures use it to determine the flat size of their variates
without inspecting each component.
"""
@inline mspace_flatsize(::Type{MU}) where {MU} = NoMSpaceElementSize{MU}()


"""
    MeasureBase.mspace_ndims(::Type{MU})
    MeasureBase.mspace_ndims(μ)

The number of dimensions of the flat variates of measures of type `MU`,
`0` for scalar variates, or a [`MeasureBase.NoMSpaceElementSize`](@ref)
if unknown.

Batched kernels rely on it to tell the variate dimensions of a flat batch
from its batch dimensions. It follows from
[`MeasureBase.mspace_flatsize`](@ref) where that is known, measure types
with array variates of dynamic size declare it directly.
"""
function mspace_ndims end

@inline mspace_ndims(::Type{MU}) where {MU} = _ndims_of_size(mspace_flatsize(MU), MU)
@inline mspace_ndims(μ::MU) where {MU} = _ndims_of_size(mspace_flatsize(μ), MU, mspace_ndims(MU))

@inline _ndims_of_size(sz::SizeLike, ::Type) = length(_size_dims(sz))
@inline _ndims_of_size(::NoMSpaceElementSize, ::Type{MU}) where {MU} = NoMSpaceElementSize{MU}()
@inline _ndims_of_size(sz::SizeLike, ::Type, ::Any) = length(_size_dims(sz))
@inline _ndims_of_size(::NoMSpaceElementSize, ::Type, n) = n

@inline _add_ndims(n::Integer, k::Integer) = n + k
@inline _add_ndims(n::NoMSpaceElementSize, ::Integer) = n
