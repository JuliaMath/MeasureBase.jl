# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsChainRulesCoreExt

using MeasureBase
import Distributions
import ChainRulesCore
using ChainRulesCore: NoTangent

using MeasureBase: _dist_params_numtype
using Distributions: Distribution

_dist_params_numtype_pullback(ΔΩ) = (NoTangent(), NoTangent())
function ChainRulesCore.rrule(::typeof(_dist_params_numtype), d::Distribution)
    _dist_params_numtype(d), _dist_params_numtype_pullback
end

end # module MeasureBaseDistributionsChainRulesCoreExt
