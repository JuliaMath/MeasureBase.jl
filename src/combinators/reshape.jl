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
struct Reshape{M<:SizeLike,N<:SizeLike} <: Function
    output_size::M
    input_size::N

    Reshape{M,N}(out_sz::M, in_sz::N) where {M<:SizeLike,N<:SizeLike} =
        new{M,N}(out_sz, in_sz)
end

function Reshape(output_size::SizeLike, input_size::SizeLike)
    out_sz = canonical_size(output_size)
    in_sz = canonical_size(input_size)
    return Reshape{typeof(out_sz), typeof(in_sz)}(out_sz, in_sz)
end

_throw_reshape_mismatch(sz, sz_x) = throw(DimensionMismatch("Reshape input size is $sz but got input of size $sz_x"))

function (f::Reshape)(x::AbstractArray)
    sz_x = maybestatic_size(x)
    f.input_size == sz_x || _throw_reshape_mismatch(f.input_size, sz_x)
    return reshape(x, f.output_size)
end

InverseFunctions.inverse(f::Reshape{M,N}) where {M,N} = Reshape{N,M}(f.input_size, f.output_size)

function ChangesOfVariables.with_logabsdet_jacobian(f::Reshape, x::AbstractArray)
    return f(x), zero(real_numtype(typeof(x)))
end


"""
    mreshape(m::AbstractMeasure, sz::Vararg{N,IntegerLike}) where N
    mreshape(m::AbstractMeasure, sz::NTuple{N,IntegerLike}) where N

Reshape a measure `m` over an array-valued space, returning a measure over
a space of arrays with shape `sz`.
"""
function mreshape end

mreshape(m::AbstractMeasure, sz::IntegerLike...) = mreshape(m, sz)
mreshape(m::AbstractMeasure, sz::SizeLike) = pushfwd(Reshape(sz, some_mspace_elsize(m)), m)

@inline mspace_elsize(μ::PushforwardMeasure{<:Reshape}) = μ.f.output_size
# Reshaped variates have flat storage only if the reshaped elements are
# numbers, the reshape of a nested variate has no flat form its density
# kernel could consume:
@inline function mspace_flatsize(μ::PushforwardMeasure{<:Reshape})
    _reshaped_flatsize(mspace_flatsize(μ.origin), mspace_elsize(μ.origin), μ.f.output_size)
end
@inline function _reshaped_flatsize(sz_flat::SizeLike, sz_outer::SizeLike, sz_out)
    _reshaped_flatsize(Val(length(_size_dims(sz_flat)) == length(_size_dims(sz_outer))), sz_out)
end
@inline _reshaped_flatsize(::Val{true}, sz_out) = sz_out
@inline _reshaped_flatsize(::Val{false}, sz_out) = NoMSpaceElementSize{typeof(sz_out)}()
@inline _reshaped_flatsize(::Any, ::Any, sz_out) = NoMSpaceElementSize{typeof(sz_out)}()
