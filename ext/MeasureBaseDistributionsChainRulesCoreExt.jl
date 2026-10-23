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
using MeasureBase: _gamma_cdf_dshape, _beta_cdf_dshape

# Derivatives of the regularized incomplete gamma and beta functions and
# their inverses: with respect to the variate resp. probability from the
# density, with respect to the shape parameters from `_gamma_cdf_dshape` and
# `_beta_cdf_dshape` (by the implicit function theorem for the quantiles):

function ChainRulesCore.rrule(::typeof(_gamma_cdf), α::Real, x::Real)
    y = _gamma_cdf(α, x)
    dy_dx = exp(_gamma_logpdf(α, x))
    function _gamma_cdf_pullback(ȳ)
        Δ = unthunk(ȳ)
        return NoTangent(), _gamma_cdf_dshape(α, x) * Δ, dy_dx * Δ
    end
    return y, _gamma_cdf_pullback
end
function ChainRulesCore.rrule(::typeof(_gamma_quantile), α::Real, p::Real)
    x = _gamma_quantile(α, p)
    dx_dp = exp(-_gamma_logpdf(α, x))
    function _gamma_quantile_pullback(x̄)
        Δ = dx_dp * unthunk(x̄)
        return NoTangent(), -_gamma_cdf_dshape(α, x) * Δ, Δ
    end
    return x, _gamma_quantile_pullback
end
function ChainRulesCore.rrule(::typeof(_beta_cdf), α::Real, β::Real, x::Real)
    y = _beta_cdf(α, β, x)
    dy_dx = exp(_beta_logpdf(α, β, x))
    function _beta_cdf_pullback(ȳ)
        Δ = unthunk(ȳ)
        dy_dα, dy_dβ = _beta_cdf_dshape(α, β, x)
        return NoTangent(), dy_dα * Δ, dy_dβ * Δ, dy_dx * Δ
    end
    return y, _beta_cdf_pullback
end
function ChainRulesCore.rrule(::typeof(_beta_quantile), α::Real, β::Real, p::Real)
    x = _beta_quantile(α, β, p)
    dx_dp = exp(-_beta_logpdf(α, β, x))
    function _beta_quantile_pullback(x̄)
        Δ = dx_dp * unthunk(x̄)
        dI_dα, dI_dβ = _beta_cdf_dshape(α, β, x)
        return NoTangent(), -dI_dα * Δ, -dI_dβ * Δ, Δ
    end
    return x, _beta_quantile_pullback
end

function ChainRulesCore.rrule(::typeof(_dist_params_numtype), d::Distribution)
    _dist_params_numtype(d), _dist_params_numtype_pullback
end

end # module MeasureBaseDistributionsChainRulesCoreExt
