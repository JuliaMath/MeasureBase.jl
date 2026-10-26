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
        ℓ = -_column_sums(abs2, Z) ./ 2 .- (_logdet_cov(d.Σ) + length(d) * log2π) / 2
        _batch_results(_masked(ℓ, _column_all(isfinite.(Xc))), X)
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

function MeasureBase.batched_transport_to_std(::Type{StdUniform}, d::Dirichlet, X::AbstractArray)
    K = length(d)
    αs, βs = _stick_breaking_params(d)
    Xc = _as_columns(X)
    # The break fractions are the shares of the components in the remaining
    # stick, given by the suffix sums (exact also for small trailing
    # components):
    S = _suffix_sums(Xc)
    z = _ratio_or_zero.(_rows(Xc, 1:(K - 1)), _rows(S, 1:(K - 1)))
    U = _beta_cdf.(αs, βs, _nonnan.(_unit_clamp.(z), 1 // 2))
    _from_columns(_nan_columns(U, _simplex_mask(Xc)), X)
end

function MeasureBase.batched_transport_from_std(::Type{StdUniform}, d::Dirichlet, U::AbstractArray)
    αs, βs = _stick_breaking_params(d)
    Uc = _as_columns(U)
    P = _quantile_arg.(Uc)
    # The break fractions `z` and their complements `1 - z`, both from the
    # lower tail like `beta_inc_inv` does, exact also for `z` close to one:
    lo = P .<= 1 // 2
    y = _beta_quantile.(ifelse.(lo, αs, βs), ifelse.(lo, βs, αs), ifelse.(lo, P, 1 .- P))
    z = ifelse.(lo, y, 1 .- y)
    r = ifelse.(lo, 1 .- y, y)
    # Each component takes its fraction of the remaining stick, the last one
    # all of it:
    X = vcat(z, _ones_row(z)) .* vcat(_ones_row(r), cumprod(r; dims = 1))
    _from_columns(_nan_columns(X, _column_all((Uc .>= 0) .& (Uc .<= 1))), U)
end

@inline _suffix_sums(X::AbstractArray) = reverse(cumsum(reverse(X; dims = 1); dims = 1); dims = 1)
@inline _ones_row(X::AbstractVector) = @_adignore fill!(similar(X, 1), 1)
@inline _ones_row(X::AbstractMatrix) = @_adignore fill!(similar(X, 1, size(X, 2)), 1)
@inline _ratio_or_zero(a, b) = ifelse(iszero(b), zero(a), a) / ifelse(iszero(b), one(b), b)

MeasureBase.transport_to_std(::Type{StdUniform}, d::Dirichlet, x) = MeasureBase.batched_transport_to_std(StdUniform, d, x)
MeasureBase.transport_from_std(::Type{StdUniform}, d::Dirichlet, z) = MeasureBase.batched_transport_from_std(StdUniform, d, z)

Adapt.adapt_structure(to, m::DirichletMeasure) = asmeasure(Dirichlet(Adapt.adapt(to, m.obj.alpha)))
