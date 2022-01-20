using Static

export PartialStatic
export partialstatic

struct PartialStatic{S,D} <: AbstractFloat
    static::S
    dynamic::D

    function PartialStatic(s::S, d::D) where {S,D}
        @assert known(is_static(s))
        new{S,D}(s, d)
    end
end

function Static.dynamic(x::PartialStatic)
    dyn = dynamic(x.dynamic)
    D = typeof(dyn)
    (x.static + dyn)::D
end

function Base.show(io::IO, p::PartialStatic)
    print(io, "PartialStatic(")
    print(io, p.static)
    print(io, ", ")
    print(io, p.dynamic)
    print(io, ")")
end

function Base.convert(::Type{PartialStatic{S, D}}, x::T) where {S, D, T<:Number}
    partialstatic(x)
end

function Base.convert(::Type{PartialStatic{S,D}}, x::T) where {S,D,T}
    partialstatic(x)
end

function Base.convert(::Type{PartialStatic{S, D}}, x::T) where {S, D, T<:(StaticFloat64)}
    partialstatic(x)
end

function Base.promote_rule(::Type{PartialStatic{S1,D1}}, ::Type{PartialStatic{S2,D2}}) where {S1,D1,S2,D2}
    PartialStatic{promote_rule(S1, S2), promote_rule(D1, D2)}
end

function promote_rule(::Type{PartialStatic{S, D}}, ::Type{StaticFloat64{N}}) where {N, S, D}
    PartialStatic{promote_rule(StaticFloat64{N}, S), D}
end

function promote_rule(::Type{StaticFloat64{N}}, ::Type{PartialStatic{S, D}}) where {N, S, D}
    PartialStatic{promote_rule(StaticFloat64{N}, S), D}
end

function Base.promote_rule(::Type{T}, ::Type{PartialStatic{S,D}}) where {T,S,D}
    _promote(PartialStatic{S,D}, T, is_static(T))
end

function Base.promote_rule(::Type{PartialStatic{S,D}}, ::Type{T}) where {S,D,T}
    _promote(PartialStatic{S,D}, T, is_static(T))
end

function _promote(::Type{PartialStatic{S,D}}, ::Type{T}, ::True) where {S,D,T}
    PartialStatic{promote_rule(S,T),D}
end

function _promote(::Type{PartialStatic{S,D}}, ::Type{T}, ::False) where {S,D,T}
    PartialStatic{S,promote_rule(D,T)}
end

partialstatic(s, d) = PartialStatic(s, d)

partialstatic(x) = _partialstatic(x, is_static(x))

partialstatic(x::PartialStatic) = x

_partialstatic(x, ::True) = PartialStatic(x, static(0.0))

_partialstatic(x, ::False) = PartialStatic(static(0.0), x)

function Base.:+(x::PartialStatic{S1,D1}, y::PartialStatic{S2,D2}) where {S1,D1,S2,D2}
    PartialStatic(x.static + y.static, x.dynamic + y.dynamic)
end

function Base.:-(x::PartialStatic{S1,D1}, y::PartialStatic{S2,D2}) where {S1,D1,S2,D2}
    PartialStatic(x.static - y.static, x.dynamic - y.dynamic)
end

Base.eps(::Type{PartialStatic{S,D}}) where {S,D} = max(Base.eps(known(S)), Base.eps(D))

# function Base.:-(x::PartialStatic{S,D}, y::StaticFloat64{Y}) where {S,D,Y}
#     PartialStatic(x.static - y, x.dynamic)
# end

# function Base.:+(x::PartialStatic{S,D}, y::Number) where {S,D}
#     PartialStatic(x.static, x.dynamic + y)
# end

# function Base.:-(x::PartialStatic{S,D}, y::Number) where {S,D}
#     PartialStatic(x.static, x.dynamic - y)
# end

# function Base.:+(x::PartialStatic{S,D}, y) where {S,D}
#     PartialStatic(x.static, x.dynamic + y)
# end

# function Base.:-(x::PartialStatic{S,D}, y) where {S,D}
#     PartialStatic(x.static, x.dynamic - y)
# end

function add_partial_static(x::StaticFloat64, y::StaticFloat64)
    return x + y
end

function add_partial_static(x::StaticFloat64, y::Real)
    return PartialStatic(x, y)
end

function add_partial_static(y::Real, x::StaticFloat64)
    return PartialStatic(x, y)
end

function add_partial_static(x, y)
    return x + y
end

Base.exp(x::PartialStatic) = exp(known(x.static) + x.dynamic)

Base.isapprox(x::PartialStatic, y::PartialStatic; kwargs...) = isapprox(dynamic(x), dynamic(y); kwargs...)
Base.isapprox(x::Number, y::PartialStatic; kwargs...) = isapprox(dynamic(x), dynamic(y); kwargs...)
Base.isapprox(x::PartialStatic, y::Number; kwargs...) = isapprox(dynamic(x), dynamic(y); kwargs...)