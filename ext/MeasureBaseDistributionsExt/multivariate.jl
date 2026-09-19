# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Multivariate normal and Dirichlet measures with density kernels and
# transports over flat batches `(n, batch dims...)`, in terms of array
# operations on the parameters, so that they run on devices and in traced
# code. Single variates are vectors.

# Column batches: flat batches as `(n, :)` matrices, single variates stay
# vectors. Reductions over the columns give arrays over the batch
# dimensions, numbers for single variates.
@inline _as_columns(x::AbstractVector) = x
@inline _as_columns(X::AbstractArray) = reshape(X, (size(X, 1), :))
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

for bhead in (:batched_logdensityof_impl, :batched_logdensity_def)
    @eval function MeasureBase.$bhead(m::MvNormalMeasure, X::AbstractArray)
        d = m.obj
        Z = _cholesky_L(d.Σ) \ (_as_columns(X) .- d.μ)
        sq = _column_sums(abs2, Z)
        _batch_results(-sq ./ 2 .- (_logdet_cov(d.Σ) + length(d) * log2π) / 2, X)
    end
end
MeasureBase.logdensity_def(m::MvNormalMeasure, x::AbstractVector) = MeasureBase.batched_logdensity_def(m, x)
MeasureBase.unsafe_logdensityof(m::MvNormalMeasure, x::AbstractVector) = MeasureBase.batched_logdensityof_impl(m, x)

function MeasureBase.batched_transport_to_std(::Type{StdNormal}, d::MvNormal, X::AbstractArray)
    _from_columns(_cholesky_L(d.Σ) \ (_as_columns(X) .- d.μ), X)
end
function MeasureBase.batched_transport_from_std(::Type{StdNormal}, d::MvNormal, Z::AbstractArray)
    _from_columns(_cholesky_L(d.Σ) * _as_columns(Z) .+ d.μ, Z)
end
MeasureBase.transport_to_std(::Type{StdNormal}, d::MvNormal, x) = MeasureBase.batched_transport_to_std(StdNormal, d, x)
MeasureBase.transport_from_std(::Type{StdNormal}, d::MvNormal, z) = MeasureBase.batched_transport_from_std(StdNormal, d, z)

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
# cumulative sums and products running along the variate dimension.

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

@inline function _simplex_mask(Xc::AbstractArray)
    tol = sqrt(eps(float(eltype(Xc))))
    _column_all(Xc .>= 0) .& (abs.(_column_sums(identity, Xc) .- 1) .<= tol)
end

# The stick-breaking Beta parameters, for the first `K - 1` components:
@inline _stick_breaking_params(d::Dirichlet) = (_dropfront(_rev_cumsum(d.alpha)), _dropback(d.alpha))

function MeasureBase.batched_transport_to_std(::Type{StdUniform}, d::Dirichlet, X::AbstractArray)
    K = length(d)
    αs, βs = _stick_breaking_params(d)
    Xc = _as_columns(X)
    rem = 1 .- cumsum(Xc; dims = 1)
    # The remaining mass before each component is the mass after it plus
    # the component itself:
    beta_v = _rows(rem, 1:(K - 1)) ./ (_rows(rem, 1:(K - 1)) .+ _rows(Xc, 1:(K - 1)))
    Z = _beta_cdf.(αs, βs, _unit_clamp.(beta_v))
    _from_columns(_nan_columns(Z, _simplex_mask(Xc)), X)
end

function MeasureBase.batched_transport_from_std(::Type{StdUniform}, d::Dirichlet, Z::AbstractArray)
    K = length(d)
    αs, βs = _stick_breaking_params(d)
    Zc = _as_columns(Z)
    beta_v = _beta_quantile.(αs, βs, _unit_clamp.(Zc))
    cp = cumprod(beta_v; dims = 1)
    # Each component takes what its Beta variate leaves of the remaining mass:
    X = vcat(1 .- _rows(cp, 1:1), _rows(cp, 1:(K - 2)) .- _rows(cp, 2:(K - 1)), _rows(cp, (K - 1):(K - 1)))
    _from_columns(_nan_columns(X, _column_all((Zc .>= 0) .& (Zc .<= 1))), Z)
end
MeasureBase.transport_to_std(::Type{StdUniform}, d::Dirichlet, x) = MeasureBase.batched_transport_to_std(StdUniform, d, x)
MeasureBase.transport_from_std(::Type{StdUniform}, d::Dirichlet, z) = MeasureBase.batched_transport_from_std(StdUniform, d, z)

Adapt.adapt_structure(to, m::DirichletMeasure) = asmeasure(Dirichlet(Adapt.adapt(to, m.obj.alpha)))
