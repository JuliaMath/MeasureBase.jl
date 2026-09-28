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

# The outermost dual number type among the arguments, duals of inner tags
# are constants with respect to it:
@inline MeasureBase._dualtag(::_Dual{T}, rest::Number...) where {T} = _outer_dual(_Dual{T}, MeasureBase._dualtag(rest...))
@inline _outer_dual(::Type{D}, ::Type{Nothing}) where {D<:_Dual} = D
@inline _outer_dual(::Type{_Dual{S}}, ::Type{_Dual{T}}) where {S,T} = S === T || !ForwardDiff.:≺(S, T) ? _Dual{S} : _Dual{T}

@inline _value(::Type{_Dual{T}}, x) where {T} = ForwardDiff.value(T, x)

# The derivative through the last argument; the partials of parameters
# are marked NaN where they are nonzero (their derivatives are unknown):
@inline _arg_partials(::Type{_Dual{T}}, x::_Dual{T}, params...) where {T} = ForwardDiff.partials(x)
@inline _arg_partials(::Type{D}, ::Real, params...) where {D<:_Dual} = zero(ForwardDiff.partials(_first_dual(D, params...)))
@inline _first_dual(::Type{_Dual{T}}, x::_Dual{T}, rest...) where {T} = x
@inline _first_dual(::Type{D}, ::Real, rest...) where {D<:_Dual} = _first_dual(D, rest...)
@inline _nan_partials(::Type{D}, ∂, ::Real) where {D<:_Dual} = ∂
@inline function _nan_partials(::Type{_Dual{T}}, ∂, x::_Dual{T}) where {T}
    ∂ + ForwardDiff.Partials(map(v -> ifelse(iszero(v), zero(v), oftype(v, NaN)), ForwardDiff.partials(x).values))
end
@inline function _through_last(::Type{_Dual{T}}, value, dvalue, last, params...) where {T}
    ∂ = dvalue * _arg_partials(_Dual{T}, last, params...)
    ForwardDiff.Dual{T}(value, foldl((∂, p) -> _nan_partials(_Dual{T}, ∂, p), params; init = ∂))
end

@inline function MeasureBase._gamma_cdf_impl(::Type{D}, α, x) where {D<:_Dual}
    αv, xv = _value(D, α), _value(D, x)
    _through_last(D, _gamma_cdf(αv, xv), exp(_gamma_logpdf(αv, xv)), x, α)
end
@inline function MeasureBase._gamma_quantile_impl(::Type{D}, α, p) where {D<:_Dual}
    αv, pv = _value(D, α), _value(D, p)
    xv = _gamma_quantile(αv, pv)
    _through_last(D, xv, exp(-_gamma_logpdf(αv, xv)), p, α)
end
@inline function MeasureBase._beta_cdf_impl(::Type{D}, α, β, x) where {D<:_Dual}
    αv, βv, xv = _value(D, α), _value(D, β), _value(D, x)
    _through_last(D, _beta_cdf(αv, βv, xv), exp(_beta_logpdf(αv, βv, xv)), x, α, β)
end
@inline function MeasureBase._beta_quantile_impl(::Type{D}, α, β, p) where {D<:_Dual}
    αv, βv, pv = _value(D, α), _value(D, β), _value(D, p)
    xv = _beta_quantile(αv, βv, pv)
    _through_last(D, xv, exp(-_beta_logpdf(αv, βv, xv)), p, α, β)
end

# The quantile of Beta doesn't support dual parameters:
@inline MeasureBase._dist_quantile(d::Beta{<:ForwardDiff.Dual}, p::Real) = convert(float(typeof(p)), NaN)
@inline MeasureBase._dist_cquantile(d::Beta{<:ForwardDiff.Dual}, p::Real) = convert(float(typeof(p)), NaN)

end # module MeasureBaseDistributionsForwardDiffExt
