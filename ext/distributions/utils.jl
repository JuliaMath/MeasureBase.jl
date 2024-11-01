# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).


"""
    convert_realtype(::Type{T}, x) where {T<:Real}

Convert x to use `T` as it's underlying type for real numbers.
"""
function convert_realtype end

_convert_realtype_pullback(ΔΩ) = NoTangent(), NoTangent, ΔΩ
ChainRulesCore.rrule(::typeof(convert_realtype), ::Type{T}, x) where T = convert_realtype(T, x), _convert_realtype_pullback

@inline convert_realtype(::Type{T}, x::T) where {T<:Real} = x
@inline convert_realtype(::Type{T}, x::AbstractArray{T}) where {T<:Real} = x
@inline convert_realtype(::Type{T}, x::U) where {T<:Real,U<:Real} = T(x)
convert_realtype(::Type{T}, x::AbstractArray{U}) where {T<:Real,U<:Real} = T.(x)
convert_realtype(::Type{T}, x) where {T<:Real} = fmap(elem -> convert_realtype(T, elem), x)


"""
    firsttype(::Type{T}, ::Type{U}) where {T<:Real,U<:Real}

Return the first type, but as a dual number type if the second one is dual.

If `U <: ForwardDiff.Dual{tag,<:Real,N}`, returns `ForwardDiff.Dual{tag,T,N}`,
otherwise returns `T`
"""
function firsttype end

firsttype(::Type{T}, ::Type{U}) where {T<:Real,U<:Real} = T
firsttype(::Type{T}, ::Type{<:ForwardDiff.Dual{tag,<:Real,N}}) where {T<:Real,tag,N} = ForwardDiff.Dual{tag,T,N}
