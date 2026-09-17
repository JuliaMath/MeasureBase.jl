# Standard measure preferences of the combinators and standard measure ranks:

@inline _stdmeasure_rank(::Type{StdUniform}) = 1
@inline _stdmeasure_rank(::Type{StdExponential}) = 2
@inline _stdmeasure_rank(::Type{StdLogistic}) = 3
@inline _stdmeasure_rank(::Type{StdNormal}) = 4


@inline preferred_stdmeasure(::Type{<:PowerMeasure{M}}) where {M} = preferred_stdmeasure(M)
@inline preferred_stdmeasure(::Type{<:WeightedMeasure{<:Any,M}}) where {M} = preferred_stdmeasure(M)
# Transports of the base don't transport the restricted measure:
@inline preferred_stdmeasure(::Type{MU}) where {MU<:RestrictedMeasure} = NoStdTransport{MU}
@inline preferred_stdmeasure(::Type{<:PushforwardMeasure{<:Any,<:Any,M}}) where {M} = preferred_stdmeasure(M)
@inline preferred_stdmeasure(::Type{<:Dirac}) = AnyStdMeasure

@inline preferred_stdmeasure(::Type{<:ProductMeasure{M}}) where {M<:AbstractArray} = preferred_stdmeasure(eltype(M))

# Arrays of marginals of mixed types combine their preferences at run time:
@inline function preferred_stdmeasure(μ::ProductMeasure{<:AbstractArray})
    _array_product_stdmeasure(μ, preferred_stdmeasure(typeof(μ)))
end
@inline _array_product_stdmeasure(μ, S::Type) = S
function _array_product_stdmeasure(μ::ProductMeasure{<:AbstractArray{M}}, ::Type{NoStdTransport{MU}}) where {M,MU}
    isconcretetype(M) && return NoStdTransport{MU}
    mapreduce(preferred_stdmeasure, promote_stdmeasure, marginals(μ); init = AnyStdMeasure)
end

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
