# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Multivariate normal and Dirichlet measures with density kernels and
# transports over flat batches `(n, batch dims...)`, in terms of array
# operations on the parameters, so that they run on devices and in traced
# code. Single variates are vectors.

# Column batches: flat batches as `(n, :)` matrices, single variates stay
# vectors. Reductions over the columns give arrays over the batch
# dimensions, numbers for single variates.
@inline _as_columns(x::AbstractVector) = x
@inline _as_columns(X::AbstractArray) = reshape(X, (size(X, 1), prod(Base.tail(size(X)))))
@inline _from_columns(y::AbstractVector, ::AbstractVector) = y
@inline _from_columns(Y::AbstractMatrix, X::AbstractArray) = reshape(Y, (size(Y, 1), Base.tail(size(X))...))
@inline _column_sums(f, z::AbstractVector) = sum(f, z)
@inline _column_sums(f, Z::AbstractMatrix) = vec(sum(f, Z; dims = 1))
@inline _column_all(z::AbstractVector) = all(z)
@inline _column_all(Z::AbstractMatrix) = vec(all(Z; dims = 1))
@inline _batch_results(r::Number, ::AbstractVector) = r
@inline _batch_results(r::AbstractVector, X::AbstractArray) = reshape(r, Base.tail(size(X)))
@inline _rows(z::AbstractVector, r) = view(z, r)
@inline _rows(Z::AbstractMatrix, r) = view(Z, r, :)
@inline _masked(ℓ::Number, ins) = ifelse(ins, ℓ, oftype(ℓ, -Inf))
@inline _masked(ℓ::AbstractArray, ins) = ifelse.(ins, ℓ, eltype(ℓ)(-Inf))
@inline _nan_columns(Z::AbstractArray, ins) = ifelse.(_as_row(ins), Z, eltype(Z)(NaN))
@inline _as_row(ins::Bool) = ins
@inline _as_row(ins::AbstractVector) = reshape(ins, 1, :)


# Multivariate normal: densities via the Cholesky factor of the covariance.

const MvNormalMeasure = AsMeasure{<:MvNormal}

_logdet_cov(Σ::PDMats.PDMat) = 2 * sum(log, diag(Σ.chol.factors))
_logdet_cov(Σ::PDMats.PDiagMat) = sum(log, Σ.diag)
_logdet_cov(Σ::PDMats.ScalMat) = Σ.dim * log(Σ.value)
_logdet_cov(Σ) = logdet(Σ)

# Infinite values lie outside the support:
for bhead in (:batched_logdensityof_impl, :batched_logdensity_def)
    @eval function MeasureBase.$bhead(m::MvNormalMeasure, X::AbstractArray)
        d = m.obj
        Xc = _as_columns(X)
        Z = _cholesky_L(d.Σ) \ (Xc .- d.μ)
        ld = _logdet_cov(d.Σ)
        ℓ = -_column_sums(abs2, Z) ./ 2 .- (ld + log2π * oftype(ld, length(d))) / 2
        _batch_results(_masked(ℓ, _column_all(isfinite.(Xc))), X)
    end
end
MeasureBase.logdensity_def(m::MvNormalMeasure, x::AbstractVector) = MeasureBase.batched_logdensity_def(m, x)
MeasureBase.unsafe_logdensityof(m::MvNormalMeasure, x::AbstractVector) = MeasureBase.batched_logdensityof_impl(m, x)

# Columns with infinite values lie outside the support:
function MeasureBase.batched_transport_to_std(::Type{StdNormal}, d::MvNormal, X::AbstractArray)
    Xc = _as_columns(X)
    Z = _cholesky_L(d.Σ) \ (Xc .- d.μ)
    _from_columns(_nan_columns(Z, _column_all(isfinite.(Xc))), X)
end
function MeasureBase.batched_transport_from_std(::Type{StdNormal}, d::MvNormal, Z::AbstractArray)
    _from_columns(_cholesky_L(d.Σ) * _as_columns(Z) .+ d.μ, Z)
end
MeasureBase.transport_to_std(::Type{StdNormal}, d::MvNormal, x) = MeasureBase.batched_transport_to_std(StdNormal, d, x)
MeasureBase.transport_from_std(::Type{StdNormal}, d::MvNormal, z) = MeasureBase.batched_transport_from_std(StdNormal, d, z)

# Finite Gaussian transports are affine. Read the current factors once per
# call, in the sample/parameter precision, rather than whitening every sample
# again in two reference densities. Other storage and AD types keep the generic
# path, including its nonfinite conventions.
const GaussianTransportMeasure = Union{MvNormalMeasure,MeasureBase.StdPowerMeasure{StdNormal}}
const GaussianTransport = MeasureBase.TransportFunction{<:GaussianTransportMeasure,<:GaussianTransportMeasure}

_gaussian_partype(m::MvNormalMeasure, T) = Distributions.partype(m.obj)
_gaussian_partype(::MeasureBase.StdPowerMeasure{StdNormal}, T) = T
_gaussian_finite(m::MvNormalMeasure) =
    m.obj.Σ isa Union{PDMats.PDMat,PDMats.PDiagMat,PDMats.ScalMat} && all(isfinite, m.obj.μ)
_gaussian_finite(::MeasureBase.StdPowerMeasure{StdNormal}) = true
_gaussian_logdet(m::MvNormalMeasure, T) = _gaussian_logdet(m.obj.Σ, T)
_gaussian_logdet(::MeasureBase.StdPowerMeasure{StdNormal}, T) = zero(T)
_gaussian_logdet(Σ::PDMats.PDMat, T) = 2 * sum(x -> log(T(x)), diag(Σ.chol.factors))
_gaussian_logdet(Σ::PDMats.PDiagMat, T) = sum(x -> log(T(x)), Σ.diag)
_gaussian_logdet(Σ::PDMats.ScalMat, T) = T(Σ.dim) * log(T(Σ.value))

function _gaussian_ladj(f::GaussianTransport, X, Y)
    Xf, Yf = MeasureBase._flat_storage(X), MeasureBase._flat_storage(Y)
    Xf isa StridedArray || return nothing
    Yf isa StridedArray || return nothing
    isconcretetype(eltype(Xf)) && eltype(Xf) <: AbstractFloat || return nothing
    T = promote_type(eltype(Xf), _gaussian_partype(f.μ, eltype(Xf)), _gaussian_partype(f.ν, eltype(Xf)))
    T <: AbstractFloat || return nothing
    _gaussian_finite(f.μ) && _gaussian_finite(f.ν) || return nothing
    all(isfinite, Xf) && all(isfinite, Yf) || return nothing
    ladj = (_gaussian_logdet(f.ν, T) - _gaussian_logdet(f.μ, T)) / 2
    return isfinite(ladj) ? ladj : nothing
end

function ChangesOfVariables.with_logabsdet_jacobian(f::GaussianTransport, x)
    y = f(x)
    ladj = _gaussian_ladj(f, x, y)
    if isnothing(ladj)
        ladj = MeasureBase._transport_ladj(logdensityof(f.μ, x), logdensityof(f.ν, y))
    end
    return y, ladj
end

function MeasureBase.transports_with_ladj(f::GaussianTransport, X)
    Y = MeasureBase.transports(f, X)
    ladj = _gaussian_ladj(f, X, Y)
    if isnothing(ladj)
        return Y, broadcast(MeasureBase._transport_ladj, MeasureBase.logdensities(f.μ, X), MeasureBase.logdensities(f.ν, Y))
    end
    dims = Base.tail(size(MeasureBase._as_stream_batch(MeasureBase._flat_batch(X), f.μ)))
    return Y, isempty(dims) ? ladj : fill(ladj, dims)
end

# Parameters follow the batches to the device:
function Adapt.adapt_structure(to, m::MvNormalMeasure)
    d = m.obj
    asmeasure(MvNormal(Adapt.adapt(to, d.μ), _adapt_cov(to, d.Σ)))
end
function _adapt_cov(to, Σ::PDMats.PDMat)
    chol = Σ.chol
    PDMats.PDMat(Adapt.adapt(to, Σ.mat), Cholesky(Adapt.adapt(to, chol.factors), chol.uplo, chol.info))
end
_adapt_cov(to, Σ::PDMats.PDiagMat) = PDMats.PDiagMat(Adapt.adapt(to, Σ.diag))
_adapt_cov(to, Σ::PDMats.ScalMat) = Σ


# Dirichlet: densities over the simplex, transports via the stick-breaking
# Beta transports (M. J. Betancourt, "Cruising The Simplex: Hamiltonian
# Monte Carlo and the Dirichlet Distribution", arXiv:1010.3436), with the
# cumulative sums and products running along the variate dimension. The
# orientation is FlatPPL's: the break fractions `z_i ~ Beta(α_i, Σ_{j>i} α_j)`
# give `x_i = z_i * (1 - Σ_{j<i} x_j)`.

for bhead in (:batched_logdensityof_impl, :batched_logdensity_def)
    @eval function MeasureBase.$bhead(m::DirichletMeasure, X::AbstractArray)
        d = m.obj
        Xc = _as_columns(X)
        ℓ = _column_sums(identity, _clog.(d.alpha .- 1, abs.(Xc))) .- d.lmnB
        _batch_results(_masked(ℓ, _simplex_mask(Xc)), X)
    end
end
MeasureBase.logdensity_def(m::DirichletMeasure, x::AbstractVector) = MeasureBase.batched_logdensity_def(m, x)
MeasureBase.unsafe_logdensityof(m::DirichletMeasure, x::AbstractVector) = MeasureBase.batched_logdensityof_impl(m, x)

@inline _simplex_mask(Xc::AbstractArray) = MeasureBase._in_stdsimplex(Xc)

# The Beta parameters of the break fractions of the first `K - 1` components:
@inline _stick_breaking_params(d::Dirichlet) = (_dropback(d.alpha), _dropfront(_rev_cumsum(d.alpha)))

# The break fractions `z` and their complements `c = 1 - z` are the shares
# of a component and of the following ones in the remaining stick, both
# from the suffix sums, so that they stay exact also for small components.
# Their incomplete beta functions give the logistic variates of both tails:
function MeasureBase.batched_transport_to_std(::Type{StdLogistic}, d::Dirichlet, X::AbstractArray)
    K = length(d)
    αs, βs = _stick_breaking_params(d)
    Xc = _as_columns(X)
    S = _suffix_sums(Xc)
    Sk = _rows(S, 1:(K - 1))
    z = _unit_clamp.(_ratio_or_zero.(_rows(Xc, 1:(K - 1)), Sk))
    c = _unit_clamp.(_ratio_or_zero.(_rows(S, 2:K), Sk))
    L = _beta_logit.(αs, βs, _nonnan.(z, 1 // 2), _nonnan.(c, 1 // 2))
    _from_columns(_nan_columns(L, _simplex_mask(Xc)), X)
end

@inline _beta_logit(α, β, z, c) = _logit_of_tails(_beta_tails(α, β, z, c)...)

function MeasureBase.batched_transport_from_std(::Type{StdLogistic}, d::Dirichlet, L::AbstractArray)
    αs, βs = _stick_breaking_params(d)
    Lc = _as_columns(L)
    # The break fractions and their complements, the smaller one of them
    # from the quantile of its tail (see `_beta_quantiles`):
    Y1 = _beta_quantile.(αs, βs, _quantile_arg.(logistic.(Lc)))
    Y2 = _beta_quantile.(βs, αs, _quantile_arg.(logistic.(.-Lc)))
    small = Y1 .<= Y2
    z = ifelse.(small, Y1, 1 .- Y2)
    r = ifelse.(small, 1 .- Y1, Y2)
    # Each component takes its fraction of the remaining stick, the last one
    # all of it:
    X = vcat(z, _ones_row(z)) .* vcat(_ones_row(r), cumprod(r; dims = 1))
    _from_columns(_nan_columns(X, _column_all(.!isnan.(Lc))), L)
end

@inline _suffix_sums(X::AbstractArray) = reverse(cumsum(reverse(X; dims = 1); dims = 1); dims = 1)
@inline _ones_row(X::AbstractVector) = @_adignore fill!(similar(X, 1), 1)
@inline _ones_row(X::AbstractMatrix) = @_adignore fill!(similar(X, 1, size(X, 2)), 1)
@inline _ratio_or_zero(a, b) = ifelse(iszero(b), zero(a), a) / ifelse(iszero(b), one(b), b)

MeasureBase.transport_to_std(::Type{StdLogistic}, d::Dirichlet, x) = MeasureBase.batched_transport_to_std(StdLogistic, d, x)
MeasureBase.transport_from_std(::Type{StdLogistic}, d::Dirichlet, l) = MeasureBase.batched_transport_from_std(StdLogistic, d, l)

Adapt.adapt_structure(to, m::DirichletMeasure) = asmeasure(Dirichlet(Adapt.adapt(to, m.obj.alpha)))
