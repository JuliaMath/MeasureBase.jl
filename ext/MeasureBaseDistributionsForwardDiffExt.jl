# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsForwardDiffExt

import MeasureBase
import Distributions
import ForwardDiff

using Distributions: Distribution, Univariate, Continuous, Beta

@inline function MeasureBase._trafo_cdf_impl(
    ::Type{<:Union{Integer,AbstractFloat}},
    d::Distribution{Univariate,Continuous},
    x::ForwardDiff.Dual{TAG},
) where {TAG}
    x_v = ForwardDiff.value(x)
    u = Distributions.cdf(d, x_v)
    dudx = Distributions.pdf(d, x_v)
    ForwardDiff.Dual{TAG}(u, dudx * ForwardDiff.partials(x))
end

@inline function MeasureBase._trafo_quantile_impl(
    ::Type{<:Union{Integer,AbstractFloat}},
    d::Distribution{Univariate,Continuous},
    u::ForwardDiff.Dual{TAG},
) where {TAG}
    x = MeasureBase._trafo_quantile_impl_generic(d, ForwardDiff.value(u))
    dxdu = inv(Distributions.pdf(d, x))
    ForwardDiff.Dual{TAG}(x, dxdu * ForwardDiff.partials(u))
end

# Workaround for Beta dist, ForwardDiff doesn't work for parameters:
@inline MeasureBase._trafo_quantile_impl_generic(d::Beta{T}, u::Real) where {T<:ForwardDiff.Dual} =
    convert(float(typeof(u)), NaN)

end # module MeasureBaseDistributionsForwardDiffExt
