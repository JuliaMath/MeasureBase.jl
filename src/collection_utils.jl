# ToDo: Add custom rrules for _split_after?

# ToDo: Use getindex instead of view for certain cases (array types)?
@inline function split_after(x::AbstractVector, n)
    i_first = firstindex(x)
    i_last = lastindex(x)
    view(x, i_first, i_first+n-1), view(x, n, i_last)
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
