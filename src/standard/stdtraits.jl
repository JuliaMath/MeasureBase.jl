"""
    struct MeasureBase.NoStdTransport{MU}

Indicates that measures of type `MU` can't be transported to or from a
standard measure.
"""
struct NoStdTransport{MU} end

"""
    struct MeasureBase.AnyStdMeasure

Indicates that any standard measure serves as transport partner, e.g.
for measures with zero degrees of freedom.
"""
struct AnyStdMeasure end

const _StdTransportPartner = Union{Type{<:StdMeasure},Type{AnyStdMeasure},Type{<:NoStdTransport}}

"""
    MeasureBase.preferred_stdmeasure(μ)::Type
    MeasureBase.preferred_stdmeasure(::Type{MU})::Type

The type of standard measure that variates of `μ` are transported to and
from by default.

Returns `MeasureBase.AnyStdMeasure` if any standard measure serves and
`MeasureBase.NoStdTransport{MU}` if measures of type `MU` have no
standard-measure transport. Composite measures combine the preferences of
their components via [`MeasureBase.promote_stdmeasure`](@ref).

Measure types that support transport to and from standard measures should
specialize the type-based method.
"""
function preferred_stdmeasure end

@inline preferred_stdmeasure(μ) = preferred_stdmeasure(typeof(μ))
@inline preferred_stdmeasure(::Type{MU}) where {MU} = NoStdTransport{MU}

@inline preferred_stdmeasure(::Type{MU}) where {MU<:StdMeasure} = MU

"""
    MeasureBase.promote_stdmeasure(A::Type, B::Type)::Type

Combine two results of [`MeasureBase.preferred_stdmeasure`](@ref) into
the preferred standard measure type of a measure composed of both.

Standard measure types promote to the one with the wider range of
values that remain distinguishable in floating point arithmetic:
`StdUniform` promotes to any other standard measure type,
`StdExponential` to `StdLogistic` and `StdNormal`, and `StdLogistic`
to `StdNormal`.
"""
function promote_stdmeasure end

@inline function promote_stdmeasure(::Type{A}, ::Type{B}) where {A<:StdMeasure,B<:StdMeasure}
    _stdmeasure_rank(A) >= _stdmeasure_rank(B) ? A : B
end

@inline promote_stdmeasure(::Type{AnyStdMeasure}, ::Type{B}) where {B} = B
@inline promote_stdmeasure(::Type{A}, ::Type{AnyStdMeasure}) where {A} = A
@inline promote_stdmeasure(::Type{AnyStdMeasure}, ::Type{AnyStdMeasure}) = AnyStdMeasure
@inline promote_stdmeasure(::Type{A}, ::Type{B}) where {A<:NoStdTransport,B} = A
@inline promote_stdmeasure(::Type{A}, ::Type{B}) where {A,B<:NoStdTransport} = B
@inline promote_stdmeasure(::Type{A}, ::Type{B}) where {A<:NoStdTransport,B<:NoStdTransport} = A
@inline promote_stdmeasure(::Type{A}, ::Type{AnyStdMeasure}) where {A<:NoStdTransport} = A
@inline promote_stdmeasure(::Type{AnyStdMeasure}, ::Type{B}) where {B<:NoStdTransport} = B

@inline promote_stdmeasure(::Type{A}) where {A} = A
@inline function promote_stdmeasure(::Type{A}, ::Type{B}, Cs::Vararg{Type,N}) where {A,B,N}
    promote_stdmeasure(promote_stdmeasure(A, B), Cs...)
end

@inline _stdmeasure_rank(::Type{<:StdMeasure}) = 0
@inline _stdmeasure_rank(::Type{StdUniform}) = 1
@inline _stdmeasure_rank(::Type{StdExponential}) = 2
@inline _stdmeasure_rank(::Type{StdLogistic}) = 3
@inline _stdmeasure_rank(::Type{StdNormal}) = 4


@inline preferred_stdmeasure(::Type{<:PowerMeasure{M}}) where {M} = preferred_stdmeasure(M)
@inline preferred_stdmeasure(::Type{<:WeightedMeasure{<:Any,M}}) where {M} = preferred_stdmeasure(M)
@inline preferred_stdmeasure(::Type{<:RestrictedMeasure{<:Any,M}}) where {M} = preferred_stdmeasure(M)
@inline preferred_stdmeasure(::Type{<:PushforwardMeasure{<:Any,<:Any,M}}) where {M} = preferred_stdmeasure(M)
@inline preferred_stdmeasure(::Type{<:Dirac}) = AnyStdMeasure

@inline preferred_stdmeasure(::Type{<:ProductMeasure{M}}) where {M<:AbstractArray} = preferred_stdmeasure(eltype(M))

@inline function preferred_stdmeasure(::Type{<:ProductMeasure{M}}) where {M<:Tuple}
    _promote_stdmeasure_oftypes(M)
end

@inline function preferred_stdmeasure(::Type{<:ProductMeasure{NamedTuple{names,M}}}) where {names,M<:Tuple}
    _promote_stdmeasure_oftypes(M)
end

@inline _promote_stdmeasure_oftypes(::Type{Tuple{}}) = AnyStdMeasure
@generated function _promote_stdmeasure_oftypes(::Type{M}) where {M<:Tuple}
    args = [:(preferred_stdmeasure($T)) for T in M.parameters]
    :(promote_stdmeasure($(args...)))
end
