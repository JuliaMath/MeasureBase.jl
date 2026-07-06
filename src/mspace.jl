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
the arrays that are the elements of the space.

May return [`NoMSpaceElementSize{typeof(μ)}()`](@ref).
"""
function mspace_elsize end
export mspace_elsize

@inline mspace_elsize(μ::AbstractMeasure) = NoMSpaceElementSize{typeof(μ)}()


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
