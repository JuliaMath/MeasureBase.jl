# ToDo: Add custom rrules for the get/view/split/etc. functions defined here.

Base.@propagate_inbounds _as_tuple(v::AbstractVector, ::Val{N}) where {N} = Tuple(SVector{N}(v))

Base.Base.@propagate_inbounds function _get_or_view(A::AbstractVector, from::IntegerLike, until::IntegerLike)
    (view(A, from:until))
end
Base.Base.@propagate_inbounds function _get_or_view(A::AbstractVector, ::StaticInteger{from}, ::StaticInteger{until}) where {from,until}
    SVector{until-from+1}(view(A, from:until))
end

# ToDo: Specialize for StaticVector instead of SVector?
Base.Base.@propagate_inbounds function _get_or_view(A::SVector, ::StaticInteger{from}, ::StaticInteger{until}) where {from,until}
    # ToDo: Improve implementation:
    SVector(_get_or_view(Tuple(A), from, until))
end

Base.Base.@propagate_inbounds function _get_or_view(tpl::Tuple, from::IntegerLike, until::IntegerLike)
    ntuple(i -> tpl[from + i - 1], Val(until - from + 1))
end
# ToDo: Is this specialization necessary?
Base.Base.@propagate_inbounds function _get_or_view(tpl::Tuple, ::StaticInteger{from}, ::StaticInteger{until}) where {from,until}
    ntuple(i -> tpl[from + i - 1], Val(until - from + 1))
end


@inline function _split_after(x::AbstractVector, n::IntegerLike)
    i_first = maybestatic_firstindex(x)
    i_last = maybestatic_lastindex(x)
    _get_or_view(x, i_first, i_first + n - static(1)), _get_or_view(x, i_first + n, i_last)
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


# ToDo: Add static reshape for static arrays!


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

# Field access functions for Fill:
_fill_value(x::Fill) = x.value
_fill_axes(x::Fill) = x.axes


_flatten_to_rv(VV::AbstractVector{<:AbstractVector{<:Real}}) = flatview(VectorOfArrays(VV))
_flatten_to_rv(VV::AbstractVector{<:StaticVector{N,<:Real}}) where N = flatview(VectorOfSimilarArrays(VV))

_flatten_to_rv(VV::VectorOfSimilarVectors{<:Real}) = flatview(VV)
_flatten_to_rv(VV::VectorOfVectors{<:Real}) = flatview(VV)

_flatten_to_rv(tpl::Tuple{<:AbstractVector{<:Real}}) = vcat(tpl...)
_flatten_to_rv(tpl::Tuple{<:StaticVector{N,<:Real}}) where N = vcat(tpl...)
