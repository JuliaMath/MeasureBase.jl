function _pushfront(v::AbstractVector, x)
    T = promote_type(eltype(v), typeof(x))
    r = similar(v, T, length(eachindex(v)) + 1)
    r[firstindex(r)] = x
    r[firstindex(r)+1:lastindex(r)] = v
    r
end

function _pushback(v::AbstractVector, x)
    T = promote_type(eltype(v), typeof(x))
    r = similar(v, T, length(eachindex(v)) + 1)
    r[lastindex(r)] = x
    r[firstindex(r):lastindex(r)-1] = v
    r
end

_dropfront(v::AbstractVector) = v[firstindex(v)+1:lastindex(v)]

_dropback(v::AbstractVector) = v[firstindex(v):lastindex(v)-1]

_rev_cumsum(xs::AbstractVector) = reverse(cumsum(reverse(xs)))

# Equivalent to `cumprod(xs)``:
_exp_cumsum_log(xs::AbstractVector) = exp.(cumsum(log.(xs)))

Base.@propagate_inbounds _as_tuple(v::AbstractVector, ::Val{N}) where {N} = Tuple(SVector{N}(v))


Base.@propagate_inbounds function _get_or_view(A::AbstractVector, from::IntegerLike, until::IntegerLike)
    view(A, dynamic(from):dynamic(until))
end

Base.@propagate_inbounds function _get_or_view(
    A::AbstractVector,
    ::StaticInteger{from},
    ::StaticInteger{until},
) where {from,until}
    SVector{until - from + 1}(view(A, from:until))
end

# ToDo: Specialize for StaticVector instead of SVector?
Base.@propagate_inbounds function _get_or_view(
    A::SVector,
    from::StaticInteger,
    until::StaticInteger,
)
    # ToDo: Improve implementation:
    SVector(_get_or_view(Tuple(A), from, until))
end

Base.@propagate_inbounds function _get_or_view(tpl::Tuple, from::IntegerLike, until::IntegerLike)
    ntuple(i -> tpl[from + i - 1], Val(until - from + 1))
end


@inline function _split_after(x::AbstractVector, n::IntegerLike)
    idxs = maybestatic_eachindex(x)
    i_first = maybestatic_first(idxs)
    i_last = maybestatic_last(idxs)
    _get_or_view(x, i_first, i_first + n - one(n)), _get_or_view(x, i_first + n, i_last)
end

@inline _split_after(x::Tuple, n) = _split_after(x::Tuple, Val{n}())
@inline _split_after(x::Tuple, ::Val{N}) where {N} = x[begin:(begin+N-1)], x[(begin+N):end]

@generated function _split_after(x::NamedTuple{names}, ::Val{names_a}) where {names,names_a}
    n = length(names_a)
    if names[begin:(begin+n-1)] == names_a
        names_b = names[(begin+n):end]
        quote
            a, b = _split_after(values(x), Val($n))
            NamedTuple{$names_a}(a), NamedTuple{$names_b}(b)
        end
    else
        quote
            throw(ArgumentError("Can't split NamedTuple{$names} after {$names_a}"))
        end
    end
end


# Field access functions for Fill:
_fill_value(x::FillArrays.Fill) = x.value
_fill_axes(x::FillArrays.Fill) = x.axes


_flatten_to_rv(VV::AbstractVector{<:AbstractVector{<:Number}}) = flatview(VectorOfArrays(VV))
_flatten_to_rv(VV::AbstractVector{<:StaticVector{N,<:Number}}) where {N} =
    flatview(VectorOfSimilarArrays(VV))

_flatten_to_rv(VV::VectorOfSimilarVectors{<:Number}) = flatview(VV)
_flatten_to_rv(VV::VectorOfVectors{<:Number}) = flatview(VV)

_flatten_to_rv(::Tuple{}) = []
_flatten_to_rv(tpl::Tuple{Vararg{AbstractVector}}) = vcat(tpl...)
_flatten_to_rv(tpl::Tuple{Vararg{StaticVector}}) = vcat(tpl...)


# Non-mutating concatenation of measure collections:
_cat_measures(a::Tuple, b::Tuple) = (a..., b...)
_cat_measures(a::AbstractVector, b::Tuple) = vcat(a, [b...])
_cat_measures(a::Tuple, b::AbstractVector) = vcat([a...], b)
_cat_measures(a::AbstractVector, b::AbstractVector) = vcat(a, b)


# Take the beginning of a flat vector stream as a variate of size `sz`,
# scalar variates have size `()` and multi-rank variates are reshaped:
Base.@propagate_inbounds _consume_from_stream(x::AbstractVector, sz::Tuple{IntegerLike}) =
    _split_after(x, sz[1])

Base.@propagate_inbounds function _consume_from_stream(x::AbstractVector, ::Tuple{})
    idxs = maybestatic_eachindex(x)
    i_first = maybestatic_first(idxs)
    x[i_first], _get_or_view(x, i_first + one(i_first), maybestatic_last(idxs))
end

function _consume_from_stream(x::AbstractVector, sz::Tuple{Vararg{IntegerLike}})
    a_flat, x_rest = _split_after(x, size2length(sz))
    return maybestatic_reshape(a_flat, sz), x_rest
end

function _consume_from_stream(x::AbstractVector, @nospecialize(sz))
    throw(ArgumentError("Can't consume a variate of size $sz from a flat vector stream"))
end
