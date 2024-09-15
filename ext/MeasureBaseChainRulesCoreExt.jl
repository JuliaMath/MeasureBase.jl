# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseChainRulesCoreExt

using MeasureBase
using ChainRulesCore: NoTangent, ZeroTangent
import ChainRulesCore


@inline function ChainRulesCore.rrule(::typeof(_checksupport), cond, result)
    y = _checksupport(cond, result)
    function _checksupport_pullback(ȳ)
        return NoTangent(), ZeroTangent(), one(ȳ)
    end
    y, _checksupport_pullback
end


_require_insupport_pullback(ΔΩ) = NoTangent(), ZeroTangent()
function ChainRulesCore.rrule(::typeof(require_insupport), μ, x)
    return require_insupport(μ, x), _require_insupport_pullback
end


_origin_depth_pullback(ΔΩ) = NoTangent(), NoTangent()
ChainRulesCore.rrule(::typeof(_origin_depth), ν) = _origin_depth(ν), _origin_depth_pullback


_check_dof_pullback(ΔΩ) = NoTangent(), NoTangent(), NoTangent()
ChainRulesCore.rrule(::typeof(check_dof), ν, μ) = check_dof(ν, μ), _check_dof_pullback


_checked_arg_pullback(ΔΩ) = NoTangent(), NoTangent(), ΔΩ
ChainRulesCore.rrule(::typeof(checked_arg), ν, x) = checked_arg(ν, x), _checked_arg_pullback


end # module MeasureBaseChainRulesCoreExt
