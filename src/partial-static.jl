using Static

export PartialStatic
export partialstatic

import Base.promote_rule

struct PartialStatic{S,D} <: AbstractFloat
    static::S
    dynamic::D

    @inline function PartialStatic(s::S, d::D) where {S,D}
        @assert known(is_static(s))
        new{S,D}(s, d)
    end
end

@inline function Static.dynamic(x::PartialStatic)
    dynamic(x.static) + dynamic(x.dynamic)
end

@inline function Base.show(io::IO, p::PartialStatic)
    print(io, "PartialStatic(")
    print(io, p.static)
    print(io, ", ")
    print(io, p.dynamic)
    print(io, ")")
end

@inline Base.convert(::Type{F}, x::PartialStatic{StaticFloat64{N}, T}) where {F<:AbstractFloat, N,T} = convert(F, N) + convert(F, x.dynamic)

@inline function Base.convert(::Type{PartialStatic{Static.StaticFloat64{N}, D}}, x::PartialStatic{Static.StaticFloat64{N}, D}) where {D, N}
    return x
end

@inline function Base.convert(::Type{PartialStatic{S, D}}, x::T) where {S, D, T<:Number}
    partialstatic(x)
end

@inline function Base.convert(::Type{PartialStatic{S,D}}, x::T) where {S,D,T}
    partialstatic(x)
end

@inline function Base.convert(::Type{PartialStatic{S, D}}, x::T) where {S, D, T<:(StaticFloat64)}
    partialstatic(x)
end

@inline function Base.promote_rule(::Type{Static.StaticFloat64{N}}, ::Type{PartialStatic{S, D}}) where {N, S, D}
    Base.promote_rule(Float64, D)
end

@inline function Base.promote_rule(::Type{PartialStatic{S1,D1}}, ::Type{PartialStatic{S2,D2}}) where {S1,D1,S2,D2}
    PartialStatic{promote_type(S1, S2), promote_rule(D1, D2)}
end

@inline function Base.promote_rule(::Type{PartialStatic{S, D}}, ::Type{StaticFloat64{N}}) where {N, S, D}
    PartialStatic{promote_type(StaticFloat64{N}, S), D}
end

@inline function Base.promote_rule(::Type{T}, ::Type{PartialStatic{S,D}}) where {T,S,D}
    _promote(PartialStatic{S,D}, T, is_static(T))
end

@inline function Base.promote_rule(::Type{PartialStatic{S,D}}, ::Type{T}) where {S,D,T}
    _promote(PartialStatic{S,D}, T, is_static(T))
end

@inline function _promote(::Type{PartialStatic{S,D}}, ::Type{T}, ::True) where {S,D,T}
    PartialStatic{promote_type(S,T),D}
end

@inline function _promote(::Type{PartialStatic{S,D}}, ::Type{T}, ::False) where {S,D,T}
    PartialStatic{S,promote_type(D,T)}
end

@inline partialstatic(s, d) = PartialStatic(s, d)

@inline partialstatic(x) = _partialstatic(x, is_static(x))

@inline partialstatic(x::PartialStatic) = x

@inline _partialstatic(x, ::True) = PartialStatic(x, static(0.0))

@inline _partialstatic(x, ::False) = PartialStatic(static(0.0), x)

@inline function Base.:+(x::Static.StaticFloat64{X}, y::PartialStatic{Static.StaticFloat64{Y}, D}) where {X,Y,D}
    PartialStatic(X+Y, y.dynamic)
end

@inline function Base.:+(x::PartialStatic{S1,D1}, y::PartialStatic{S2,D2}) where {S1,D1,S2,D2}
    PartialStatic(x.static + y.static, x.dynamic + y.dynamic)
end

@inline function Base.:-(x::PartialStatic{S1,D1}, y::PartialStatic{S2,D2}) where {S1,D1,S2,D2}
    PartialStatic(x.static - y.static, x.dynamic - y.dynamic)
end

Base.eps(::Type{PartialStatic{S,D}}) where {S,D} = max(Base.eps(known(S)), Base.eps(D))

# @inline function Base.:-(x::PartialStatic{S,D}, y::StaticFloat64{Y}) where {S,D,Y}
#     PartialStatic(x.static - y, x.dynamic)
# end

# @inline function Base.:+(x::PartialStatic{S,D}, y::Number) where {S,D}
#     PartialStatic(x.static, x.dynamic + y)
# end

# @inline function Base.:-(x::PartialStatic{S,D}, y::Number) where {S,D}
#     PartialStatic(x.static, x.dynamic - y)
# end

# @inline function Base.:+(x::PartialStatic{S,D}, y) where {S,D}
#     PartialStatic(x.static, x.dynamic + y)
# end

# @inline function Base.:-(x::PartialStatic{S,D}, y) where {S,D}
#     PartialStatic(x.static, x.dynamic - y)
# end

@inline function add_partial_static(x::StaticFloat64, y::StaticFloat64)
    return x + y
end

@inline function add_partial_static(x::StaticFloat64, y::Real)
    return PartialStatic(x, y)
end

@inline function add_partial_static(y::Real, x::StaticFloat64)
    return PartialStatic(x, y)
end

@inline function add_partial_static(x, y)
    return x + y
end

Base.exp(x::PartialStatic) = exp(known(x.static) + x.dynamic)

@inline Base.isapprox(x::PartialStatic, y::PartialStatic; kwargs...) = isapprox(dynamic(x), dynamic(y); kwargs...)
@inline Base.isapprox(x::Number, y::PartialStatic; kwargs...) = isapprox(dynamic(x), dynamic(y); kwargs...)
@inline Base.isapprox(x::PartialStatic, y::Number; kwargs...) = isapprox(dynamic(x), dynamic(y); kwargs...)