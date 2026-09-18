# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Densities and standard transports of the main distribution families,
# implemented as plain arithmetic on the parameters without branches or
# foreign calls, so that the batched kernels of the wrapped measures run on
# devices and in traced code. Distributions' own implementations remain in
# use for other families.

const _Families = Union{Normal,Uniform,Exponential,Logistic,Cauchy,Laplace,LogNormal,Weibull,Gamma,Beta,Poisson,Bernoulli}

# Densities relative to the base measures (Lebesgue resp. counting
# measure), support checks are separate masks. The density formulas must
# not throw outside of the support, where their results are masked:
@inline MeasureBase.logdensity_def(m::AsMeasure{<:_Families}, x) = _family_logd(m.obj, x)
@inline MeasureBase.unsafe_logdensityof(m::AsMeasure{<:_Families}, x) = _family_logd(m.obj, x)
@inline MeasureBase.insupport(m::AsMeasure{<:_Families}, x) = _family_insupport(m.obj, x)

# `c * log(y)`, zero for `c == 0` also where `y == 0`:
@inline _clog(c, y) = ifelse(iszero(c), zero(c * log(one(y))), c * log(y))

@inline function _family_logd(d::Normal, x)
    z = (x - d.μ) / d.σ
    -z * z / 2 - log(d.σ) - log2π / 2
end
@inline _family_insupport(::Normal, x) = true

@inline _family_logd(d::Uniform, x) = -log(d.b - d.a) + zero(x)
@inline _family_insupport(d::Uniform, x) = (d.a <= x) & (x <= d.b)

@inline _family_logd(d::Exponential, x) = -x / d.θ - log(d.θ)
@inline _family_insupport(::Exponential, x) = x >= 0

@inline function _family_logd(d::Logistic, x)
    z = (x - d.μ) / d.θ
    -z - 2 * log1pexp(-z) - log(d.θ)
end
@inline _family_insupport(::Logistic, x) = true

@inline function _family_logd(d::Cauchy, x)
    z = (x - d.μ) / d.σ
    -log1p(z * z) - log(π * d.σ)
end
@inline _family_insupport(::Cauchy, x) = true

@inline _family_logd(d::Laplace, x) = -abs((x - d.μ) / d.θ) - log(2 * d.θ)
@inline _family_insupport(::Laplace, x) = true

@inline function _family_logd(d::LogNormal, x)
    lx = log(abs(x))
    z = (lx - d.μ) / d.σ
    ℓ = -z * z / 2 - log(d.σ) - log2π / 2 - lx
    ifelse(x > 0, ℓ, oftype(ℓ, -Inf))
end
@inline _family_insupport(::LogNormal, x) = x >= 0

@inline function _family_logd(d::Weibull, x)
    xθ = abs(x / d.θ)
    ℓ = log(d.α / d.θ) + _clog(d.α - 1, xθ) - xθ^d.α
    ifelse(isinf(xθ), oftype(ℓ, -Inf), ℓ)
end
@inline _family_insupport(::Weibull, x) = x >= 0

@inline function _family_logd(d::Gamma, x)
    ℓ = _clog(d.α - 1, abs(x)) - x / d.θ - loggamma(d.α) - d.α * log(d.θ)
    ifelse(isinf(x), oftype(ℓ, -Inf), ℓ)
end
@inline _family_insupport(::Gamma, x) = x >= 0

@inline function _family_logd(d::Beta, x)
    _clog(d.α - 1, abs(x)) + _clog(d.β - 1, abs(1 - x)) - logbeta(d.α, d.β)
end
@inline _family_insupport(::Beta, x) = (0 <= x) & (x <= 1)

@inline _family_logd(d::Poisson, x) = _clog(x, d.λ) - d.λ - loggamma(abs(x) + 1)
@inline _family_insupport(::Poisson, x) = (x >= 0) & (x == floor(x))

@inline _family_logd(d::Bernoulli, x) = ifelse(x == 1, log(d.p), log1p(-d.p))
@inline _family_insupport(::Bernoulli, x) = (x == 0) | (x == 1)


# Standard transports of the non-affine families: Cauchy, Laplace, Gamma
# and Beta pivot on the uniform measure, the log-normal and Weibull
# families on the normal resp. exponential measure.

# The regularized incomplete gamma and beta functions and their inverses
# (from SpecialFunctions), with derivatives with respect to the variate
# resp. probability argument provided by the autodiff extensions:
@inline MeasureBase._gamma_cdf(α, x) = MeasureBase._gamma_cdf_impl(MeasureBase._dualtag(α, x), α, x)
@inline MeasureBase._gamma_quantile(α, p) = MeasureBase._gamma_quantile_impl(MeasureBase._dualtag(α, p), α, p)
@inline MeasureBase._beta_cdf(α, β, x) = MeasureBase._beta_cdf_impl(MeasureBase._dualtag(α, β, x), α, β, x)
@inline MeasureBase._beta_quantile(α, β, p) = MeasureBase._beta_quantile_impl(MeasureBase._dualtag(α, β, p), α, β, p)
@inline MeasureBase._gamma_cdf_impl(::Type{Nothing}, α, x) = first(gamma_inc(α, x))
# The complementary probability is formed in the common float type, as
# `gamma_inc_inv` requires `p + q == 1` exactly:
@inline function MeasureBase._gamma_quantile_impl(::Type{Nothing}, α, p)
    T = float(promote_type(typeof(α), typeof(p)))
    pp = convert(T, p)
    gamma_inc_inv(convert(T, α), pp, one(T) - pp)
end
@inline MeasureBase._beta_cdf_impl(::Type{Nothing}, α, β, x) = first(beta_inc(α, β, x))
@inline MeasureBase._beta_quantile_impl(::Type{Nothing}, α, β, p) = first(beta_inc_inv(α, β, p))
@inline MeasureBase._gamma_logpdf(α, x) = _clog(α - 1, x) - x - loggamma(α)
@inline MeasureBase._beta_logpdf(α, β, x) = _clog(α - 1, x) + _clog(β - 1, 1 - x) - logbeta(α, β)

@inline MeasureBase.preferred_stdmeasure(::Type{<:Cauchy}) = StdUniform
@inline MeasureBase.transport_to_std(::Type{StdUniform}, d::Cauchy, x) = 1 // 2 + atan((x - d.μ) / d.σ) / π
@inline MeasureBase.transport_from_std(::Type{StdUniform}, d::Cauchy, p) = muladd(d.σ, tan(π * (p - 1 // 2)), d.μ)

@inline MeasureBase.preferred_stdmeasure(::Type{<:Laplace}) = StdUniform
@inline function MeasureBase.transport_to_std(::Type{StdUniform}, d::Laplace, x)
    z = (x - d.μ) / d.θ
    ifelse(z < 0, exp(z) / 2, 1 - exp(-z) / 2)
end
@inline function MeasureBase.transport_from_std(::Type{StdUniform}, d::Laplace, p)
    u = p - 1 // 2
    muladd(-d.θ * sign(u), log1p(-2 * abs(u)), d.μ)
end

@inline MeasureBase.preferred_stdmeasure(::Type{<:LogNormal}) = StdNormal
@inline MeasureBase.transport_to_std(::Type{StdNormal}, d::LogNormal, x) = (log(x) - d.μ) / d.σ
@inline MeasureBase.transport_from_std(::Type{StdNormal}, d::LogNormal, z) = exp(muladd(d.σ, z, d.μ))

@inline MeasureBase.preferred_stdmeasure(::Type{<:Weibull}) = StdExponential
@inline MeasureBase.transport_to_std(::Type{StdExponential}, d::Weibull, x) = (x / d.θ)^d.α
@inline MeasureBase.transport_from_std(::Type{StdExponential}, d::Weibull, z) = d.θ * z^(1 / d.α)

@inline MeasureBase.preferred_stdmeasure(::Type{<:Gamma}) = StdUniform
@inline MeasureBase.transport_to_std(::Type{StdUniform}, d::Gamma, x) = _gamma_cdf(d.α, x / d.θ)
@inline MeasureBase.transport_from_std(::Type{StdUniform}, d::Gamma, p) = d.θ * _gamma_quantile(d.α, p)

@inline MeasureBase.preferred_stdmeasure(::Type{<:Beta}) = StdUniform
@inline MeasureBase.transport_to_std(::Type{StdUniform}, d::Beta, x) = _beta_cdf(d.α, d.β, x)
@inline MeasureBase.transport_from_std(::Type{StdUniform}, d::Beta, p) = _beta_quantile(d.α, d.β, p)
