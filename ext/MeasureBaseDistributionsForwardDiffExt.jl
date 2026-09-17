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

# The quantile of Beta doesn't support dual parameters:
@inline MeasureBase._dist_quantile(d::Beta{<:ForwardDiff.Dual}, p::Real) = convert(float(typeof(p)), NaN)
@inline MeasureBase._dist_cquantile(d::Beta{<:ForwardDiff.Dual}, p::Real) = convert(float(typeof(p)), NaN)

end # module MeasureBaseDistributionsForwardDiffExt
