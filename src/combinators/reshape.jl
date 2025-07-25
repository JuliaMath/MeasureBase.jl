# ToDo: Support static resizes for static arrays

"""
    struct MeasureBase.Reshape <: Function

Represents a function that reshapes an array.

Supports `InverseFunctions.inverse` and
`ChangesOfVariables.with_logabsdet_jacobian`.

Constructor:

```julia
Reshape(output_size::Dims, input_size::Dims)
```
"""
struct Reshape{M,N} <: Function
    output_size::NTuple{M,Int}
    input_size::NTuple{N,Int}
end

_throw_reshape_mismatch(sz, sz_x) = throw(DimensionMismatch("Reshape input size is $sz but got input of size $sz_x"))

function (f::Reshape)(x::AbstractArray)
    sz_x = size(x)
    f.input_size == sz_x || _throw_reshape_mismatch(f.input_size, sz_x)
    return reshape(x, f.output_size)
end

InverseFunctions.inverse(f::Reshape) = Reshape(f.input_size, f.output_size)

ChangesOfVariables.with_logabsdet_jacobian(::Reshape, x::AbstractArray) = zero(real_numtype(typeof(x)))


"""
    mreshape(m::AbstractMeasure, sz::Vararg{N,Integer}) where N
    mreshape(m::AbstractMeasure, sz::NTuple{N,Integer}) where N

Reshape a measure `m` over an array-valued space, returning a measure over
a space of arrays with shape `sz`.
"""
function mreshape end

_elsize_for_reshape(m::AbstractMeasure) = _elsize_for_reshape(some_mspace_elsize(m), m)
_elsize_for_reshape(sz::NTuple{<:Any,Integer}, ::AbstractMeasure) = sz
#!!!!!!!!!!_elsize_for_reshape(::NoMSpaceElementSize, m::AbstractMeasure) = size(testvalue(m))

mreshape(m::AbstractMeasure, sz::Vararg{Any,N}) where N = mreshape(m, sz)
mreshape(m::AbstractMeasure, sz::NTuple{N,<:Any}) where N = pushfwd(Reshape(sz, _elsize_for_reshape(m)), m)
