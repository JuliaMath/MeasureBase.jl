# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseReactantExt

using Reactant: TracedRNumber
import MeasureBase
using MeasureBase: RealValues, IntegerValues

Base.in(::TracedRNumber{<:Real}, ::RealValues) = true
Base.in(::TracedRNumber{<:Integer}, ::IntegerValues) = true

end # module MeasureBaseReactantExt
