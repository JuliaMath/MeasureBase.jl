# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsEnzymeCoreExt

using EnzymeCore: EnzymeRules
using MeasureBase: _gamma_cdf, _gamma_ccdf, _gamma_quantile, _gamma_cquantile
using MeasureBase: _beta_cdf, _beta_quantile_pair, _gamma_logpdf, _beta_logpdf
using MeasureBase: _gamma_cdf_dshape, _beta_cdf_dshape

# Differentiate the scalar kernels analytically, including the inverse
# functions by implicit differentiation, rather than their iterative solvers.
EnzymeRules.@easy_rule(
    _gamma_cdf(α::Real, x::Real),
    (_gamma_cdf_dshape(α, x), exp(_gamma_logpdf(α, x))),
)

EnzymeRules.@easy_rule(
    _gamma_ccdf(α::Real, x::Real),
    (-_gamma_cdf_dshape(α, x), -exp(_gamma_logpdf(α, x))),
)

EnzymeRules.@easy_rule(
    _gamma_quantile(α::Real, p::Real),
    @setup(dx_dp = exp(-_gamma_logpdf(α, Ω))),
    (-_gamma_cdf_dshape(α, Ω) * dx_dp, dx_dp),
)

EnzymeRules.@easy_rule(
    _gamma_cquantile(α::Real, q::Real),
    @setup(dx_dq = -exp(-_gamma_logpdf(α, Ω))),
    (_gamma_cdf_dshape(α, Ω) * dx_dq, dx_dq),
)

# Compute the joint shape partials only when needed.
EnzymeRules.@easy_rule(
    _beta_cdf(α::Real, β::Real, x::Real),
    @setup(dshape = nothing),
    (
        first(@something(dshape, (dshape = _beta_cdf_dshape(α, β, x)))),
        last(@something(dshape, (dshape = _beta_cdf_dshape(α, β, x)))),
        exp(_beta_logpdf(α, β, x)),
    ),
)

EnzymeRules.@easy_rule(
    _beta_quantile_pair(α::Real, β::Real, p::Real),
    @setup((x, y) = Ω, dshape = nothing, dx_dp = exp(-_beta_logpdf(α, β, x, y))),
    (
        -first(@something(dshape, (dshape = _beta_cdf_dshape(α, β, x, y)))) * dx_dp,
        -last(@something(dshape, (dshape = _beta_cdf_dshape(α, β, x, y)))) * dx_dp,
        dx_dp,
    ),
    (
        first(@something(dshape, (dshape = _beta_cdf_dshape(α, β, x, y)))) * dx_dp,
        last(@something(dshape, (dshape = _beta_cdf_dshape(α, β, x, y)))) * dx_dp,
        -dx_dp,
    ),
)

end # module MeasureBaseDistributionsEnzymeCoreExt
