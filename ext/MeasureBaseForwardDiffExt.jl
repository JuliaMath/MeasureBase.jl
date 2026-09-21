# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseForwardDiffExt

using MeasureBase
using MeasureBase: containsnan
import ForwardDiff

function MeasureBase.containsnan(x::ForwardDiff.Dual)
    a = containsnan(x.value)
    b = containsnan(x.partials)
    return a || b
end


end # module MeasureBaseForwardDiffExt
