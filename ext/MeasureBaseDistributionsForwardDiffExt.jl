# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsForwardDiffExt

import MeasureBase
import Distributions
import ForwardDiff

using Distributions: Distribution, Univariate, Continuous, Beta

# Dual-number transports for distributions with plain parameters, via the
# derivatives of cdf and quantile:

const _PlainParams = Type{<:Union{Integer,AbstractFloat}}

@inline function MeasureBase._trafo_logcdf_impl(
    ::_PlainParams,
    d::Distribution{Univariate,Continuous},
    x::ForwardDiff.Dual{TAG},
) where {TAG}
    x_v = ForwardDiff.value(x)
    lp = Distributions.logcdf(d, x_v)
    dlp_dx = exp(Distributions.logpdf(d, x_v) - lp)
    ForwardDiff.Dual{TAG}(lp, dlp_dx * ForwardDiff.partials(x))
end

@inline function MeasureBase._trafo_logccdf_impl(
    ::_PlainParams,
    d::Distribution{Univariate,Continuous},
    x::ForwardDiff.Dual{TAG},
) where {TAG}
    x_v = ForwardDiff.value(x)
    lp = Distributions.logccdf(d, x_v)
    dlp_dx = -exp(Distributions.logpdf(d, x_v) - lp)
    ForwardDiff.Dual{TAG}(lp, dlp_dx * ForwardDiff.partials(x))
end

@inline function MeasureBase._trafo_quantile_impl(
    ::_PlainParams,
    d::Distribution{Univariate,Continuous},
    p::ForwardDiff.Dual{TAG},
) where {TAG}
    p_v = ForwardDiff.value(p)
    x = MeasureBase._dist_quantile(d, p_v)
    dx_dp = inv(Distributions.pdf(d, x))
    ForwardDiff.Dual{TAG}(x, dx_dp * ForwardDiff.partials(p))
end

@inline function MeasureBase._trafo_cquantile_impl(
    ::_PlainParams,
    d::Distribution{Univariate,Continuous},
    p::ForwardDiff.Dual{TAG},
) where {TAG}
    p_v = ForwardDiff.value(p)
    x = MeasureBase._dist_cquantile(d, p_v)
    dx_dp = -inv(Distributions.pdf(d, x))
    ForwardDiff.Dual{TAG}(x, dx_dp * ForwardDiff.partials(p))
end

# Dual numbers through the regularized incomplete gamma and beta functions
# and their inverses: the derivative with respect to the variate resp.
# probability argument follows from the density, derivatives with respect
# to the parameters are not available and yield NaN partials.
using MeasureBase: _gamma_cdf, _gamma_quantile, _beta_cdf, _beta_quantile, _gamma_logpdf, _beta_logpdf

const _Dual = ForwardDiff.Dual

@inline MeasureBase._dualtag(::_Dual{TAG}, ::Number...) where {TAG} = _Dual{TAG}

# The derivative through the last argument; the partials of parameters
# are marked NaN where they are nonzero (their derivatives are unknown):
@inline _arg_partials(::Type{_Dual{TAG}}, x::_Dual{TAG}, params...) where {TAG} = ForwardDiff.partials(x)
@inline _arg_partials(::Type{_Dual{TAG}}, ::Real, params...) where {TAG} = zero(ForwardDiff.partials(_first_dual(params...)))
@inline _first_dual(x::_Dual, rest...) = x
@inline _first_dual(::Real, rest...) = _first_dual(rest...)
@inline _nan_partials(∂, ::Real) = ∂
@inline function _nan_partials(∂, x::_Dual)
    ∂ + ForwardDiff.Partials(map(v -> ifelse(iszero(v), zero(v), oftype(v, NaN)), ForwardDiff.partials(x).values))
end
@inline function _through_last(::Type{_Dual{TAG}}, value, dvalue, last, params...) where {TAG}
    ∂ = dvalue * _arg_partials(_Dual{TAG}, last, params...)
    ForwardDiff.Dual{TAG}(value, foldl(_nan_partials, params; init = ∂))
end

@inline function MeasureBase._gamma_cdf_impl(::Type{_Dual{TAG}}, α, x) where {TAG}
    αv, xv = ForwardDiff.value(α), ForwardDiff.value(x)
    _through_last(_Dual{TAG}, _gamma_cdf(αv, xv), exp(_gamma_logpdf(αv, xv)), x, α)
end
@inline function MeasureBase._gamma_quantile_impl(::Type{_Dual{TAG}}, α, p) where {TAG}
    αv, pv = ForwardDiff.value(α), ForwardDiff.value(p)
    xv = _gamma_quantile(αv, pv)
    _through_last(_Dual{TAG}, xv, exp(-_gamma_logpdf(αv, xv)), p, α)
end
@inline function MeasureBase._beta_cdf_impl(::Type{_Dual{TAG}}, α, β, x) where {TAG}
    αv, βv, xv = ForwardDiff.value(α), ForwardDiff.value(β), ForwardDiff.value(x)
    _through_last(_Dual{TAG}, _beta_cdf(αv, βv, xv), exp(_beta_logpdf(αv, βv, xv)), x, α, β)
end
@inline function MeasureBase._beta_quantile_impl(::Type{_Dual{TAG}}, α, β, p) where {TAG}
    αv, βv, pv = ForwardDiff.value(α), ForwardDiff.value(β), ForwardDiff.value(p)
    xv = _beta_quantile(αv, βv, pv)
    _through_last(_Dual{TAG}, xv, exp(-_beta_logpdf(αv, βv, xv)), p, α, β)
end

# The quantile of Beta doesn't support dual parameters:
@inline MeasureBase._dist_quantile(d::Beta{<:ForwardDiff.Dual}, p::Real) = convert(float(typeof(p)), NaN)
@inline MeasureBase._dist_cquantile(d::Beta{<:ForwardDiff.Dual}, p::Real) = convert(float(typeof(p)), NaN)

end # module MeasureBaseDistributionsForwardDiffExt
