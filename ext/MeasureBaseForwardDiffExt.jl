# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseForwardDiffExt

using MeasureBase
using MeasureBase: containsnan, firsttype
import ForwardDiff

function MeasureBase.containsnan(x::ForwardDiff.Dual)
    a = containsnan(x.value)
    b = containsnan(x.partials)
    return a || b
end

MeasureBase.firsttype(::Type{T}, ::Type{<:ForwardDiff.Dual{tag,<:Real,N}}) where {T<:Real,tag,N} =
    ForwardDiff.Dual{tag,T,N}

end # module MeasureBaseForwardDiffExt
