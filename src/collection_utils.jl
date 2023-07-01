# ToDo: Specialize for StaticArray:
@inline _getindex_or_view(A, idxs...) = view(A, idxs...)


# ToDo: Add custom rrules for _split_after?

# ToDo: Use getindex instead of view for certain cases (array types)?
@inline function _split_after(x::AbstractVector, n)
    i_first = firstindex(x)
    i_last = lastindex(x)
    _getindex_or_view(x, i_first, i_first+n-1), _getindex_or_view(x, n, i_last)
end

@inline _split_after(x::Tuple, n) = _split_after(x::Tuple, Val{n}())
@inline _split_after(x::Tuple, ::Val{N}) where N = x[begin:begin+N-1], x[N:end]

@generated function _split_after(x::NamedTuple{names}, ::Val{names_a}) where {names, names_a}
    n = length(names_a)
    if names_after[begin:begin+n-1] == names_a
        names_b = names[n:end]
        quote
            a, b = _split_after(x, Val(n))
            NamedTuple{names_a}(a), NamedTuple{names_b}(b)
        end
    else
        quote
            throw(ArgumentError("Can't split NamedTuple{$names} after {$names_a}"))
        end
    end
end


Base.@propagate_inbounds function _as_tuple(v::AbstractVector, ::Val{N}) where {N}
    @boundcheck @assert length(v) == N # ToDo: Throw proper exception
    ntuple(i -> v[i], Val(N))
end


_empty_zero(::AbstractVector{T}) where {T<:Real} = Fill(zero(T), 0)


struct _TupleNamer{names} <: Function end
(::TupleNamer{names})(x::Tuple) where names = NamedTuple{names}(x)
InverseFunctions.inverse(::TupleNamer{names}) where names = TupleUnNamer{names}()
ChangesOfVariables.with_logabsdet_jacobian(::TupleNamer{names}, x::Tuple) where names = static(false)

struct _TupleUnNamer{names} <: Function end
(::TupleUnNamer{names})(x::NamedTuple{names}) where {names} = values(x)
InverseFunctions.inverse(::TupleUnNamer{names}) where names = TupleNamer{names}()
ChangesOfVariables.with_logabsdet_jacobian(::TupleUnNamer{names}, x::NamedTuple{names}) where names = static(false)


_reorder_nt(x::NamedTuple{names},::Val{names}) where {names} = x

@generated function _reorder_nt(x::NamedTuple{names},::Val{new_names}) where {names,new_names}
    if sort([names...]) != sort([new_names...])
        :(throw(ArgumentError("Can't reorder NamedTuple{$names} to NamedTuple{$new_names}")))
    else
        expr = :(())
        for nm in new_names
            push!(expr.args, :($nm = x.$nm))
        end
        return expr
    end
end

# ToDo: Add custom rrule for _reorder_nt?


_fill_value(x::Fill) = x.value
_fill_axes(x::Fill) = x.axes


# ToDo: Flatten to SVector:
_flatten_to_rv(tpl::Tuple{Vararg{Number}}) = vcat(x...)

# ToDo:
#_flatten_to_vector(tpl::Tuple{Vararg{StaticVector}}) = ...

_flatten_to_flatten_to_numvectorvector(tpl::Tuple{AbstractVector{<:Number}}) = vcat(tpl...)

# ToDo: Use faster implementation that handles large numbers of vectors efficiently:
_flatten_to_rv(V::AbstractVector{<:AbstractVector{<:Number}}) = vcat(V...)

_flatten_to_rv(V::ArrayOfSimilarArray{<:Number}) = flatview(A)
