# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsForwardDiffPullbacksExt

import MeasureBase
using MeasureBase: StdMeasure, transport_to_std, transport_from_std

import Distributions
using Distributions: Distribution, Univariate

import ChainRulesCore
using ForwardDiffPullbacks: fwddiff

# Use ForwardDiff for univariate transports:
@inline function ChainRulesCore.rrule(::typeof(transport_to_std), ::Type{S}, d::Distribution{Univariate}, x::Any) where {S<:StdMeasure}
    ChainRulesCore.rrule(fwddiff(transport_to_std), S, d, x)
end
@inline function ChainRulesCore.rrule(::typeof(transport_from_std), ::Type{S}, d::Distribution{Univariate}, z::Any) where {S<:StdMeasure}
    ChainRulesCore.rrule(fwddiff(transport_from_std), S, d, z)
end

end # module MeasureBaseDistributionsForwardDiffPullbacksExt
