# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsForwardDiffPullbacksExt

import MeasureBase
using MeasureBase: StdMeasure, transport_def

import Distributions
using Distributions: Distribution, Univariate

import ChainRulesCore
using ForwardDiffPullbacks: fwddiff

# Use ForwardDiff for univariate transformations:
@inline function ChainRulesCore.rrule(::typeof(transport_def), ν::Distribution{Univariate}, μ::Distribution{Univariate}, x::Any)
    ChainRulesCore.rrule(fwddiff(transport_def), ν, μ, x)
end
@inline function ChainRulesCore.rrule(::typeof(transport_def), ν::StdMeasure, μ::Distribution{Univariate}, x::Any)
    ChainRulesCore.rrule(fwddiff(transport_def), ν, μ, x)
end
@inline function ChainRulesCore.rrule(::typeof(transport_def), ν::Distribution{Univariate}, μ::StdMeasure, x::Any)
    ChainRulesCore.rrule(fwddiff(transport_def), ν, μ, x)
end

end # module MeasureBaseDistributionsForwardDiffPullbacksExt
