# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsChainRulesCoreExt

using MeasureBase
import Distributions
import ChainRulesCore
using ChainRulesCore: NoTangent, ProjectTo, unthunk

using MeasureBase: _dist_params_numtype
using Distributions: Distribution

_dist_params_numtype_pullback(ΔΩ) = (NoTangent(), NoTangent())
using MeasureBase: _gamma_cdf, _gamma_ccdf, _gamma_quantile, _gamma_cquantile, _beta_cdf, _gamma_logpdf, _beta_logpdf
using MeasureBase: _gamma_cdf_dshape, _beta_cdf_dshape, _beta_quantile_pair

# Derivatives of the regularized incomplete gamma and beta functions and
# their inverses: with respect to the variate resp. probability from the
# density, with respect to the shape parameters from `_gamma_cdf_dshape` and
# `_beta_cdf_dshape` (by the implicit function theorem for the quantiles).
# The cotangents are projected onto the argument types, which may differ in
# precision:

function ChainRulesCore.rrule(::typeof(_gamma_cdf), α::Real, x::Real)
    y = _gamma_cdf(α, x)
    dy_dx = exp(_gamma_logpdf(α, x))
    function _gamma_cdf_pullback(ȳ)
        Δ = unthunk(ȳ)
        return NoTangent(), ProjectTo(α)(_gamma_cdf_dshape(α, x) * Δ), ProjectTo(x)(dy_dx * Δ)
    end
    return y, _gamma_cdf_pullback
end
function ChainRulesCore.rrule(::typeof(_gamma_quantile), α::Real, p::Real)
    x = _gamma_quantile(α, p)
    dx_dp = exp(-_gamma_logpdf(α, x))
    function _gamma_quantile_pullback(x̄)
        Δ = dx_dp * unthunk(x̄)
        return NoTangent(), ProjectTo(α)(-_gamma_cdf_dshape(α, x) * Δ), ProjectTo(p)(Δ)
    end
    return x, _gamma_quantile_pullback
end
function ChainRulesCore.rrule(::typeof(_gamma_ccdf), α::Real, x::Real)
    y = _gamma_ccdf(α, x)
    dy_dx = -exp(_gamma_logpdf(α, x))
    function _gamma_ccdf_pullback(ȳ)
        Δ = unthunk(ȳ)
        return NoTangent(), ProjectTo(α)(-_gamma_cdf_dshape(α, x) * Δ), ProjectTo(x)(dy_dx * Δ)
    end
    return y, _gamma_ccdf_pullback
end
function ChainRulesCore.rrule(::typeof(_gamma_cquantile), α::Real, q::Real)
    x = _gamma_cquantile(α, q)
    dx_dq = -exp(-_gamma_logpdf(α, x))
    function _gamma_cquantile_pullback(x̄)
        Δ = dx_dq * unthunk(x̄)
        return NoTangent(), ProjectTo(α)(_gamma_cdf_dshape(α, x) * Δ), ProjectTo(q)(Δ)
    end
    return x, _gamma_cquantile_pullback
end
function ChainRulesCore.rrule(::typeof(_beta_cdf), α::Real, β::Real, x::Real)
    y = _beta_cdf(α, β, x)
    dy_dx = exp(_beta_logpdf(α, β, x))
    function _beta_cdf_pullback(ȳ)
        Δ = unthunk(ȳ)
        dy_dα, dy_dβ = _beta_cdf_dshape(α, β, x)
        return NoTangent(), ProjectTo(α)(dy_dα * Δ), ProjectTo(β)(dy_dβ * Δ), ProjectTo(x)(dy_dx * Δ)
    end
    return y, _beta_cdf_pullback
end
function ChainRulesCore.rrule(::typeof(_beta_quantile_pair), α::Real, β::Real, p::Real)
    x, y = _beta_quantile_pair(α, β, p)
    dx_dp = exp(-_beta_logpdf(α, β, x, y))
    function _beta_quantile_pullback(pair̄)
        d = unthunk(pair̄)
        Δ = if d isa ChainRulesCore.AbstractZero
            d
        else
            unthunk(first(d)) - unthunk(last(d))
        end
        Δ = dx_dp * Δ
        dI_dα, dI_dβ = _beta_cdf_dshape(α, β, x, y)
        return (
            NoTangent(),
            ProjectTo(α)(-dI_dα * Δ),
            ProjectTo(β)(-dI_dβ * Δ),
            ProjectTo(p)(Δ),
        )
    end
    return (x, y), _beta_quantile_pullback
end

function ChainRulesCore.rrule(::typeof(_dist_params_numtype), d::Distribution)
    _dist_params_numtype(d), _dist_params_numtype_pullback
end

end # module MeasureBaseDistributionsChainRulesCoreExt
