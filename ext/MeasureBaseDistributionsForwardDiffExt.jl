# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsForwardDiffExt

import MeasureBase
import Distributions
import ForwardDiff

using Distributions: Distribution, Univariate, Continuous, Beta

# Dual-number transports for distributions with plain parameters, via the
# derivatives of cdf and quantile. The values come from the hooks again, so
# that nested dual numbers peel off one layer at a time:

const _PlainParams = Type{<:Union{Integer,AbstractFloat}}

@inline MeasureBase._positive_finite(x::ForwardDiff.Dual) =
    MeasureBase._positive_finite(ForwardDiff.value(x))

@inline function MeasureBase._trafo_logcdf_impl(
    T::_PlainParams,
    d::Distribution{Univariate,Continuous},
    x::ForwardDiff.Dual{TAG},
) where {TAG}
    x_v = ForwardDiff.value(x)
    lp = MeasureBase._trafo_logcdf_impl(T, d, x_v)
    dlp_dx = exp(Distributions.logpdf(d, x_v) - lp)
    ForwardDiff.Dual{TAG}(lp, dlp_dx * ForwardDiff.partials(x))
end

@inline function MeasureBase._trafo_logccdf_impl(
    T::_PlainParams,
    d::Distribution{Univariate,Continuous},
    x::ForwardDiff.Dual{TAG},
) where {TAG}
    x_v = ForwardDiff.value(x)
    lp = MeasureBase._trafo_logccdf_impl(T, d, x_v)
    dlp_dx = -exp(Distributions.logpdf(d, x_v) - lp)
    ForwardDiff.Dual{TAG}(lp, dlp_dx * ForwardDiff.partials(x))
end

@inline function MeasureBase._trafo_quantile_impl(
    T::_PlainParams,
    d::Distribution{Univariate,Continuous},
    p::ForwardDiff.Dual{TAG},
) where {TAG}
    p_v = ForwardDiff.value(p)
    x = MeasureBase._trafo_quantile_impl(T, d, p_v)
    dx_dp = inv(Distributions.pdf(d, x))
    ForwardDiff.Dual{TAG}(x, dx_dp * ForwardDiff.partials(p))
end

@inline function MeasureBase._trafo_cquantile_impl(
    T::_PlainParams,
    d::Distribution{Univariate,Continuous},
    p::ForwardDiff.Dual{TAG},
) where {TAG}
    p_v = ForwardDiff.value(p)
    x = MeasureBase._trafo_cquantile_impl(T, d, p_v)
    dx_dp = -inv(Distributions.pdf(d, x))
    ForwardDiff.Dual{TAG}(x, dx_dp * ForwardDiff.partials(p))
end

# Stieltjes measure functions and their inverses of distributions with plain
# parameters at dual numbers, which Distributions doesn't support for all
# families (e.g. TDist), via the hooks above:
@inline function MeasureBase._dist_call(f::F, d::Distribution{Univariate,Continuous}, x::ForwardDiff.Dual) where {F}
    _dual_call(f, MeasureBase._dist_params_numtype(d), d, x)
end
@inline _dual_call(f::F, ::Type, d, x) where {F} = f(d, x)
@inline _dual_call(::typeof(Distributions.logcdf), T::_PlainParams, d, x) = MeasureBase._trafo_logcdf_impl(T, d, x)
@inline _dual_call(::typeof(Distributions.logccdf), T::_PlainParams, d, x) = MeasureBase._trafo_logccdf_impl(T, d, x)
@inline _dual_call(::typeof(Distributions.cdf), T::_PlainParams, d, x) = exp(MeasureBase._trafo_logcdf_impl(T, d, x))
@inline _dual_call(::typeof(Distributions.ccdf), T::_PlainParams, d, x) = exp(MeasureBase._trafo_logccdf_impl(T, d, x))
@inline _dual_call(::typeof(MeasureBase._dist_quantile), T::_PlainParams, d, p) = MeasureBase._trafo_quantile_impl(T, d, p)
@inline _dual_call(::typeof(MeasureBase._dist_cquantile), T::_PlainParams, d, p) = MeasureBase._trafo_cquantile_impl(T, d, p)

# Dual numbers through the regularized incomplete gamma and beta functions
# and their inverses: the derivatives with respect to the variate resp.
# probability follow from the density, the ones with respect to the shape
# parameters from `_gamma_cdf_dshape` and `_beta_cdf_dshape` (by the
# implicit function theorem for the quantiles). The latter are only
# computed for shape parameters with partials of the tag.
using MeasureBase: _gamma_cdf, _gamma_ccdf, _gamma_quantile, _gamma_cquantile, _beta_cdf, _gamma_logpdf, _beta_logpdf
using MeasureBase: _gamma_cdf_dshape, _beta_cdf_dshape, _beta_quantile_pair

const _Dual = ForwardDiff.Dual

# The outermost dual number type among the arguments, duals of inner tags
# are constants with respect to it:
@inline MeasureBase._dualtag(::_Dual{T}, rest::Number...) where {T} = _outer_dual(_Dual{T}, MeasureBase._dualtag(rest...))
@inline _outer_dual(::Type{D}, ::Type{Nothing}) where {D<:_Dual} = D
@inline _outer_dual(::Type{_Dual{S}}, ::Type{_Dual{T}}) where {S,T} = S === T || !ForwardDiff.:≺(S, T) ? _Dual{S} : _Dual{T}

@inline _value(::Type{_Dual{T}}, x) where {T} = ForwardDiff.value(T, x)

@inline _isdual(::Type{D}, xs...) where {D<:_Dual} = any(x -> x isa D, xs)

# The dual number with the given value and the partials ∑ᵢ ∂ᵢ ẋᵢ over the
# arguments xᵢ with partials ẋᵢ of the tag:
@inline function _chain(::Type{_Dual{T}}, value, ∂s::Tuple, xs::Tuple) where {T}
    ForwardDiff.Dual{T}(value, reduce(_add_partials, map((∂, x) -> _scaled_partials(_Dual{T}, ∂, x), ∂s, xs)))
end
@inline _scaled_partials(::Type{_Dual{T}}, ∂, x::_Dual{T}) where {T} = ∂ * ForwardDiff.partials(x)
@inline _scaled_partials(::Type{<:_Dual}, ∂, ::Real) = nothing
@inline _add_partials(a, b) = a + b
@inline _add_partials(a, ::Nothing) = a
@inline _add_partials(::Nothing, b) = b
@inline _add_partials(::Nothing, ::Nothing) = nothing

@inline function MeasureBase._gamma_cdf_impl(::Type{D}, α, x) where {D<:_Dual}
    αv, xv = _value(D, α), _value(D, x)
    ∂α = _isdual(D, α) ? _gamma_cdf_dshape(αv, xv) : nothing
    _chain(D, _gamma_cdf(αv, xv), (∂α, exp(_gamma_logpdf(αv, xv))), (α, x))
end
@inline function MeasureBase._gamma_quantile_impl(::Type{D}, α, p) where {D<:_Dual}
    αv, pv = _value(D, α), _value(D, p)
    xv = _gamma_quantile(αv, pv)
    ∂p = exp(-_gamma_logpdf(αv, xv))
    ∂α = _isdual(D, α) ? -_gamma_cdf_dshape(αv, xv) * ∂p : nothing
    _chain(D, xv, (∂α, ∂p), (α, p))
end
@inline function MeasureBase._gamma_ccdf_impl(::Type{D}, α, x) where {D<:_Dual}
    αv, xv = _value(D, α), _value(D, x)
    ∂α = _isdual(D, α) ? -_gamma_cdf_dshape(αv, xv) : nothing
    _chain(D, _gamma_ccdf(αv, xv), (∂α, -exp(_gamma_logpdf(αv, xv))), (α, x))
end
@inline function MeasureBase._gamma_cquantile_impl(::Type{D}, α, q) where {D<:_Dual}
    αv, qv = _value(D, α), _value(D, q)
    xv = _gamma_cquantile(αv, qv)
    ∂q = -exp(-_gamma_logpdf(αv, xv))
    ∂α = _isdual(D, α) ? _gamma_cdf_dshape(αv, xv) * ∂q : nothing
    _chain(D, xv, (∂α, ∂q), (α, q))
end
@inline function MeasureBase._beta_cdf_impl(::Type{D}, α, β, x) where {D<:_Dual}
    αv, βv, xv = _value(D, α), _value(D, β), _value(D, x)
    ∂α, ∂β = _isdual(D, α, β) ? _beta_cdf_dshape(αv, βv, xv) : (nothing, nothing)
    _chain(D, _beta_cdf(αv, βv, xv), (∂α, ∂β, exp(_beta_logpdf(αv, βv, xv))), (α, β, x))
end
@inline function MeasureBase._beta_quantile_pair_impl(::Type{D}, α, β, p) where {D<:_Dual}
    αv, βv, pv = _value(D, α), _value(D, β), _value(D, p)
    x, y = _beta_quantile_pair(αv, βv, pv)
    ∂p = exp(-_beta_logpdf(αv, βv, x, y))
    ∂α, ∂β = if _isdual(D, α, β)
        map(∂ -> -∂ * ∂p, _beta_cdf_dshape(αv, βv, x, y))
    else
        (nothing, nothing)
    end
    xd = _chain(D, x, (∂α, ∂β, ∂p), (α, β, p))
    (xd, D(y, -ForwardDiff.partials(xd)))
end

# The quantile of Beta doesn't support dual parameters:
@inline MeasureBase._dist_quantile(d::Beta{<:ForwardDiff.Dual}, p::Real) = convert(float(typeof(p)), NaN)
@inline MeasureBase._dist_cquantile(d::Beta{<:ForwardDiff.Dual}, p::Real) = convert(float(typeof(p)), NaN)

end # module MeasureBaseDistributionsForwardDiffExt
