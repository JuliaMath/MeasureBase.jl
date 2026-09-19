# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseAffineMapsExt

using MeasureBase
using MeasureBase: PushforwardMeasure, AdaptRootMeasure, PushfwdRootMeasure
using MeasureBase: StaticInteger
using AffineMaps: AbstractAffineMap
using ChangesOfVariables: with_logabsdet_jacobian

# Affine maps treat matrices as batches of column vectors, so flat batches
# of vector variates are applied as `(n, :)` matrices:
function MeasureBase._apply_generic(f::AbstractAffineMap, X::AbstractArray, ::StaticInteger{1})
    _columns_back(f(_as_columns(X)), X)
end

@inline _as_columns(x::AbstractVector) = x
@inline _as_columns(X::AbstractArray) = reshape(X, (size(X, 1), :))
@inline _columns_back(y::AbstractVector, ::AbstractVector) = y
@inline _columns_back(Y::AbstractMatrix, X::AbstractArray) = reshape(Y, (size(Y, 1), Base.tail(size(X))...))

const _AffinePushfwd{M,S} = PushforwardMeasure{<:AbstractAffineMap,<:AbstractAffineMap,M,S}

# Densities of affine pushforwards of vector variates use the per-column
# log-abs-det-Jacobians of the inverse map:
for (bhead, head) in [(:batched_logdensityof_impl, :logdensityof_impl), (:batched_logdensity_def, :logdensity_def)]
    @eval function MeasureBase.$bhead(ν::_AffinePushfwd{M,<:AdaptRootMeasure}, Y::AbstractArray) where {M}
        _affine_pushfwd_ld(MeasureBase.$head, ν, Y, MeasureBase._static_ndims(ν))
    end
    @eval function MeasureBase.$bhead(ν::_AffinePushfwd{M,<:PushfwdRootMeasure}, Y::AbstractArray) where {M}
        MeasureBase._batched_kernel(MeasureBase.$head, ν.origin, MeasureBase._apply_batched(ν.finv, Y, MeasureBase._static_ndims(ν)))
    end
end

function _affine_pushfwd_ld(f::F, ν::PushforwardMeasure, Y::AbstractArray, ::StaticInteger{1}) where {F}
    X2, ladj2 = with_logabsdet_jacobian(ν.finv, _as_columns(Y))
    ℓ = MeasureBase._batched_kernel(f, ν.origin, _columns_back(X2, Y))
    return MeasureBase._lazy_combine_ladj(ℓ, _ladj_back(ladj2, Y))
end
_affine_pushfwd_ld(f::F, ν::PushforwardMeasure, Y::AbstractArray, k) where {F} = MeasureBase._default_batched_kernel(f, ν, Y, k)

@inline _ladj_back(ladj::Number, ::AbstractVector) = ladj
@inline _ladj_back(ladj::AbstractMatrix, Y::AbstractArray) = reshape(ladj, Base.tail(size(Y)))

end # module MeasureBaseAffineMapsExt
