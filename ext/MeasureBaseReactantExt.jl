# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseReactantExt

using Reactant: Reactant, TracedRNumber
using IrrationalConstants: sqrt2
import MeasureBase
using MeasureBase: RealValues, IntegerValues

Base.in(::TracedRNumber{<:Real}, ::RealValues) = true
Base.in(::TracedRNumber{<:Integer}, ::IntegerValues) = true

# CHLO provides erf_inv but no erfc_inv, so the standard normal quantile
# loses precision for arguments close to 0 and 1 in traced code:
MeasureBase.Φinv(p::TracedRNumber) = Reactant.Ops.erf_inv(2 * p - 1) * sqrt2

end # module MeasureBaseReactantExt
