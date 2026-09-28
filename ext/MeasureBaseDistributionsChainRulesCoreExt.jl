# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsChainRulesCoreExt

using MeasureBase
import Distributions
import ChainRulesCore
using ChainRulesCore: NoTangent, unthunk

using MeasureBase: _dist_params_numtype
using Distributions: Distribution

_dist_params_numtype_pullback(ΔΩ) = (NoTangent(), NoTangent())
using MeasureBase: _gamma_cdf, _gamma_quantile, _beta_cdf, _beta_quantile, _gamma_logpdf, _beta_logpdf

# Derivatives with respect to the variate resp. probability argument of the
# regularized incomplete gamma and beta functions and their inverses. The
# derivatives with respect to the shape parameters are not implemented:
const _NoShapeDerivative = "derivatives with respect to shape parameters of incomplete gamma and beta functions"
_∂shape() = ChainRulesCore.@not_implemented(_NoShapeDerivative)

function ChainRulesCore.rrule(::typeof(_gamma_cdf), α::Real, x::Real)
    y = _gamma_cdf(α, x)
    dy_dx = exp(_gamma_logpdf(α, x))
    return y, ȳ -> (NoTangent(), _∂shape(), dy_dx * unthunk(ȳ))
end
function ChainRulesCore.rrule(::typeof(_gamma_quantile), α::Real, p::Real)
    x = _gamma_quantile(α, p)
    dx_dp = exp(-_gamma_logpdf(α, x))
    return x, x̄ -> (NoTangent(), _∂shape(), dx_dp * unthunk(x̄))
end
function ChainRulesCore.rrule(::typeof(_beta_cdf), α::Real, β::Real, x::Real)
    y = _beta_cdf(α, β, x)
    dy_dx = exp(_beta_logpdf(α, β, x))
    return y, ȳ -> (NoTangent(), _∂shape(), _∂shape(), dy_dx * unthunk(ȳ))
end
function ChainRulesCore.rrule(::typeof(_beta_quantile), α::Real, β::Real, p::Real)
    x = _beta_quantile(α, β, p)
    dx_dp = exp(-_beta_logpdf(α, β, x))
    return x, x̄ -> (NoTangent(), _∂shape(), _∂shape(), dx_dp * unthunk(x̄))
end

function ChainRulesCore.rrule(::typeof(_dist_params_numtype), d::Distribution)
    _dist_params_numtype(d), _dist_params_numtype_pullback
end

end # module MeasureBaseDistributionsChainRulesCoreExt
