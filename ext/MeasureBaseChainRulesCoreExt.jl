# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseChainRulesCoreExt

using MeasureBase
using ChainRulesCore: NoTangent, ZeroTangent
import ChainRulesCore

# = utils ====================================================================

using MeasureBase: isneginf, isposinf

_isneginf_pullback(::Any) = (NoTangent(), ZeroTangent())
ChainRulesCore.rrule(::typeof(isneginf), x) = isneginf(x), _logdensityof_rt_pullback

_isposinf_pullback(::Any) = (NoTangent(), ZeroTangent())
ChainRulesCore.rrule(::typeof(isposinf), x) = isposinf(x), _isposinf_pullback

# = collection utils =========================================================

using MeasureBase: _dropfront, _dropback, _rev_cumsum, _exp_cumsum_log

function ChainRulesCore.rrule(::typeof(_pushfront), v::AbstractVector, x)
    result = _pushfront(v, x)
    function _pushfront_pullback(thunked_ΔΩ)
        ΔΩ = ChainRulesCore.unthunk(thunked_ΔΩ)
        (NoTangent(), ΔΩ[firstindex(ΔΩ)+1:lastindex(ΔΩ)], ΔΩ[firstindex(ΔΩ)])
    end
    return result, _pushfront_pullback
end


function ChainRulesCore.rrule(::typeof(_pushback), v::AbstractVector, x)
    result = _pushback(v, x)
    function _pushback_pullback(thunked_ΔΩ)
        ΔΩ = ChainRulesCore.unthunk(thunked_ΔΩ)
        (NoTangent(), ΔΩ[firstindex(ΔΩ):lastindex(ΔΩ)-1], ΔΩ[lastindex(ΔΩ)])
    end
    return result, _pushback_pullback
end


function ChainRulesCore.rrule(::typeof(_rev_cumsum), xs::AbstractVector)
    result = _rev_cumsum(xs)
    function _rev_cumsum_pullback(ΔΩ)
        ∂xs = ChainRulesCore.@thunk cumsum(ChainRulesCore.unthunk(ΔΩ))
        (NoTangent(), ∂xs)
    end
    return result, _rev_cumsum_pullback
end


function ChainRulesCore.rrule(::typeof(_exp_cumsum_log), xs::AbstractVector)
    result = _exp_cumsum_log(xs)
    function _exp_cumsum_log_pullback(ΔΩ)
        ∂xs = inv.(xs) .* _rev_cumsum(exp.(cumsum(log.(xs))) .* ChainRulesCore.unthunk(ΔΩ))
        (NoTangent(), ∂xs)
    end
    return result, _exp_cumsum_log_pullback
end


# = insupport & friends ======================================================

using MeasureBase: check_dof, require_insupport, checked_arg, _checksupport, _origin_depth

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

# = return type inference ====================================================

using MeasureBase: logdensityof_rt

_logdensityof_rt_pullback(::Any) = (NoTangent(), NoTangent(), ZeroTangent())
function ChainRulesCore.rrule(::typeof(logdensityof_rt), target, v)
    logdensityof_rt(target, v), _logdensityof_rt_pullback
end

end # module MeasureBaseChainRulesCoreExt
