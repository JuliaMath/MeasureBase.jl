# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseAffineMapsExt

using MeasureBase
using MeasureBase: PushforwardMeasure, AdaptRootMeasure, PushfwdRootMeasure
using MeasureBase: StaticInteger
using AffineMaps: AbstractAffineMap, Mul, Add, Subtract, MulAdd, AddMul, InvMul, InvMulAdd, InvAddMul
using ChangesOfVariables: with_logabsdet_jacobian
using LinearAlgebra: diag
import Statistics

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
# log-abs-det-Jacobians of the inverse map. Pushforwards that adapt the root
# measure have the root measure of their origin as their base measure, so
# both kernels use the full density of the origin:
for bhead in (:batched_logdensityof_impl, :batched_logdensity_def)
    @eval function MeasureBase.$bhead(ν::_AffinePushfwd{M,<:AdaptRootMeasure}, Y::AbstractArray) where {M}
        _affine_pushfwd_ld(ν, Y, MeasureBase._stream_ndims(ν))
    end
end
for (bhead, head) in [(:batched_logdensityof_impl, :logdensityof_impl), (:batched_logdensity_def, :logdensity_def)]
    @eval function MeasureBase.$bhead(ν::_AffinePushfwd{M,<:PushfwdRootMeasure}, Y::AbstractArray) where {M}
        MeasureBase._batched_kernel(MeasureBase.$head, ν.origin, MeasureBase._apply_batched(ν.finv, Y, MeasureBase._stream_ndims(ν)))
    end
end

function _affine_pushfwd_ld(ν::PushforwardMeasure, Y::AbstractArray, ::StaticInteger{1})
    X2, ladj2 = with_logabsdet_jacobian(ν.finv, _as_columns(Y))
    ℓ = MeasureBase.batched_logdensityof_impl(ν.origin, _columns_back(X2, Y))
    return MeasureBase._lazy_combine_ladj(ℓ, _ladj_back(ladj2, Y))
end
_affine_pushfwd_ld(ν::PushforwardMeasure, Y::AbstractArray, k) = MeasureBase._default_batched_kernel(MeasureBase.logdensityof_impl, ν, Y, k)

@inline _ladj_back(ladj::Number, ::AbstractVector) = ladj
@inline _ladj_back(ladj::AbstractMatrix, Y::AbstractArray) = reshape(ladj, Base.tail(size(Y)))


# Moments of affine pushforwards follow from those of their origin,
# (co)variances for scalar and vector variates:
Statistics.mean(ν::_AffinePushfwd) = ν.f(Statistics.mean(ν.origin))
Statistics.var(ν::_AffinePushfwd) = _affine_var(ν, Statistics.var(ν.origin))
Statistics.cov(ν::_AffinePushfwd) = _affine_cov(ν, Statistics.var(ν.origin))

_affine_var(ν, v::Number) = _cov_through(ν.f, v)
_affine_var(ν, ::AbstractVector{<:Number}) = diag(Statistics.cov(ν))
_affine_var(ν, ::Any) = MeasureBase._no_moment(ν, :variance)
_affine_cov(ν, ::AbstractVector{<:Number}) = _cov_through(ν.f, Statistics.cov(ν.origin))
_affine_cov(ν, ::Any) = MeasureBase._no_moment(ν, :covariance)

_cov_through(f::Union{Mul,MulAdd,AddMul}, C) = f.A * C * f.A'
_cov_through(::Union{Add,Subtract}, C) = C
_cov_through(f::Union{InvMul,InvMulAdd,InvAddMul}, C) = (f.A \ C) / f.A'

end # module MeasureBaseAffineMapsExt
