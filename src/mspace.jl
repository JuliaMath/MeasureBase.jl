"""
    MeasureBase.NoMSpaceElementSize{MU}

Indicates that there either the measurable space of measures of type `MU`
is not a space over arrays, of that the size of the arrays is not fixed
or that it can not be easily/efficiently determined.
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

@inline mspace_elsize(μ::AbstraceMeasure) = NoMSpaceElementSize{typeof(μ)}()


"""
    some_mspace_elsize(μ::AbstractMeasure)

For a measure `μ` over an array-valued measurable space, return the size of
an arbitrary element of the space.

Use with caution, the space of some measures is made up of arrays of different
sizes!

In general, use [`mspace_elsize(μ)`](@ref) instead. `some_mspace_elsize` is
useful if you the measurable space is expected to contains only arrays of the
same size but if there is no way to prove this automatically. Algorithms that
uses the returned size should always check that it matches the size of each
points of the space that is processed.
"""
function some_mspace_elsize end
export some_mspace_elsize

@inline some_mspace_elsize(μ) = _mspace_some_elsize_impl(μ, mspace_elsize(μ))

@inline _mspace_some_elsize_impl(::AbstractMeasure, sz::NTuple{<:Any,Integer}) = sz
_mspace_some_elsize_impl(μ::AbstractMeasure, ::NoMSpaceElementSize) = size(testvalue(μ))
