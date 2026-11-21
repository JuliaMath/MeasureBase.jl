# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsForwardDiffPullbacksExt

import MeasureBase
using MeasureBase: StdLogistic, _to_logistic_via, _from_logistic_via

using Distributions: Distribution, Univariate, Continuous

import ChainRulesCore
using ChainRulesCore: NoTangent, Tangent, @thunk, unthunk
using ForwardDiffPullbacks: fwddiff

# The generic transports via Distributions' cdfs and quantiles (the families
# with branch-free kernels differentiate natively) use ForwardDiff. The
# cotangent of the distribution is computed lazily per field, so that
# variate derivatives also work for distributions built inside the
# differentiated function from constant parameters ForwardDiff can't
# differentiate:
for f in (:_to_logistic_via, :_from_logistic_via)
    @eval function ChainRulesCore.rrule(::typeof($f), ::Type{StdLogistic}, d::Distribution{Univariate,Continuous}, x::Real)
        y, fwddiff_pullback = ChainRulesCore.rrule(fwddiff($f), StdLogistic, d, x)
        function $(Symbol(f, :_pullback))(Δy)
            _, _, Δd, Δx = fwddiff_pullback(Δy)
            return NoTangent(), NoTangent(), _lazy_field_tangents(d, Δd), Δx
        end
        return y, $(Symbol(f, :_pullback))
    end
end

function _lazy_field_tangents(::D, Δd) where {D}
    names = fieldnames(D)
    Tangent{D}(; NamedTuple{names}(map(n -> @thunk(getproperty(unthunk(Δd), n)), names))...)
end

end # module MeasureBaseDistributionsForwardDiffPullbacksExt
