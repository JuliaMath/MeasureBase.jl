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
@inline MeasureBase.insupportof(m::AsMeasure{<:_Families}, x) = _family_insupport(m.obj, x) & _finite_variate(m.obj, x)

# Zero by value: dual numbers compare lexicographically with their
# partials, but not across the smallest normal float (subnormals count as
# zero, as on devices that flush them):
@inline _iszero_value(x) = abs(x) < floatmin(float(typeof(x)))

# `c * log(y)` for `y >= 0`, zero for `c == y == 0`. Only that endpoint is
# guarded, so that the derivative with respect to `c` survives at `c == 0`:
@inline _clog(c, y) = c * log(ifelse(_iszero_value(c) & _iszero_value(y), one(y), y))

@inline function _family_logd(d::Normal, x)
    z = (x - d.μ) / d.σ
    -z * z / 2 - log(d.σ) - log2π / 2
end
@inline _family_insupport(::Normal, x) = true

@inline _family_logd(d::Uniform, x) = -log(d.b - d.a) + zero(x)
@inline _family_insupport(d::Uniform, x) = (d.a <= x) & (x <= d.b)

@inline _family_logd(d::Exponential, x) = -x / d.θ - log(d.θ)
@inline _family_insupport(::Exponential, x) = x >= 0

# The symmetric form of the standard logistic density stays finite in both tails:
@inline _family_logd(d::Logistic, x) = MeasureBase.logdensityof_impl(StdLogistic(), (x - d.μ) / d.θ) - log(d.θ)
@inline _family_insupport(::Logistic, x) = true

# `log1p(z^2)` as `2 log(a) + log1p((b / a)^2)` for `a = max(|z|, 1)` and
# `b = min(|z|, 1)`, so that `z^2` can't overflow:
@inline function _family_logd(d::Cauchy, x)
    z = abs((x - d.μ) / d.σ)
    big = z > 1
    a, b = ifelse(big, z, one(z)), ifelse(big, one(z), z)
    -2 * log(a) - log1p(abs2(b / a)) - log(π * d.σ)
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
# (from SpecialFunctions), differentiable with respect to all arguments via
# the autodiff extensions.
@inline MeasureBase._gamma_cdf(α, x) = MeasureBase._gamma_cdf_impl(MeasureBase._dualtag(α, x), α, x)
@inline MeasureBase._gamma_ccdf(α, x) = MeasureBase._gamma_ccdf_impl(MeasureBase._dualtag(α, x), α, x)
@inline MeasureBase._gamma_quantile(α, p) = MeasureBase._gamma_quantile_impl(MeasureBase._dualtag(α, p), α, p)
@inline MeasureBase._gamma_cquantile(α, q) = MeasureBase._gamma_cquantile_impl(MeasureBase._dualtag(α, q), α, q)
@inline MeasureBase._beta_cdf(α, β, x) = MeasureBase._beta_cdf_impl(MeasureBase._dualtag(α, β, x), α, β, x)
@inline MeasureBase._beta_quantile(α, β, p) = MeasureBase._beta_quantile_impl(MeasureBase._dualtag(α, β, p), α, β, p)
@inline MeasureBase._gamma_cdf_impl(::Type{Nothing}, α, x) = first(gamma_inc(α, x))
@inline MeasureBase._gamma_ccdf_impl(::Type{Nothing}, α, x) = last(gamma_inc(α, x))
# The complementary probability is formed in the common float type, as
# `gamma_inc_inv` requires `p + q == 1` exactly:
@inline function MeasureBase._gamma_quantile_impl(::Type{Nothing}, α, p)
    T = float(promote_type(typeof(α), typeof(p)))
    pp = convert(T, p)
    gamma_inc_inv(convert(T, α), pp, one(T) - pp)
end
@inline function MeasureBase._gamma_cquantile_impl(::Type{Nothing}, α, q)
    T = float(promote_type(typeof(α), typeof(q)))
    qq = convert(T, q)
    gamma_inc_inv(convert(T, α), one(T) - qq, qq)
end
@inline MeasureBase._beta_cdf_impl(::Type{Nothing}, α, β, x) = first(beta_inc(α, β, x))
@inline MeasureBase._beta_quantile_impl(::Type{Nothing}, α, β, p) = first(beta_inc_inv(α, β, p))
@inline MeasureBase._gamma_logpdf(α, x) = _clog(α - 1, x) - x - loggamma(α)
@inline MeasureBase._beta_logpdf(α, β, x) = _clog(α - 1, x) + _clog(β - 1, 1 - x) - logbeta(α, β)

# Branch-free formulas evaluate all branches, also in AD pullbacks, so
# unselected branches must stay finite with finite derivatives. Tails are
# computed from tail probabilities, selected by the same condition as the
# result, so that derivatives at the median follow one branch.

@inline MeasureBase.preferred_stdmeasure(::Type{<:Cauchy}) = StdUniform
# `atan(z)/π + 1/2` cancels in the lower tail, the angle `atan(1, -z)` doesn't:
@inline function MeasureBase.transport_to_std(::Type{StdUniform}, d::Cauchy, x)
    z = (x - d.μ) / d.σ
    _nan_outside(d, x, atan(one(z), -z) / π)
end
# `tan(π(p - 1/2))` loses small probabilities in `p - 1/2`, the cotangent
# of the tail probability doesn't:
@inline function MeasureBase.transport_from_std(::Type{StdUniform}, d::Cauchy, p)
    q = _unit_interior(p)
    lower = q < 1 // 2
    r = ifelse(lower, q, 1 - q)
    t = cospi(r) / sinpi(r)
    _nan_outside(StdUniform, p, muladd(d.σ, ifelse(lower, -t, t), d.μ))
end

@inline MeasureBase.preferred_stdmeasure(::Type{<:Laplace}) = StdUniform
@inline function MeasureBase.transport_to_std(::Type{StdUniform}, d::Laplace, x)
    z = (x - d.μ) / d.θ
    lower = z < 0
    t = exp(ifelse(lower, z, -z)) / 2
    _nan_outside(d, x, ifelse(lower, t, 1 - t))
end
@inline function MeasureBase.transport_from_std(::Type{StdUniform}, d::Laplace, p)
    q = _unit_interior(p)
    lower = q < 1 // 2
    l = log(2 * ifelse(lower, q, 1 - q))
    _nan_outside(StdUniform, p, muladd(d.θ, ifelse(lower, l, -l), d.μ))
end

@inline MeasureBase.preferred_stdmeasure(::Type{<:LogNormal}) = StdNormal
@inline MeasureBase.transport_to_std(::Type{StdNormal}, d::LogNormal, x) = _nan_outside(d, x, (log(abs(x)) - d.μ) / d.σ)
@inline MeasureBase.transport_from_std(::Type{StdNormal}, d::LogNormal, z) = exp(muladd(d.σ, z, d.μ))

@inline MeasureBase.preferred_stdmeasure(::Type{<:Weibull}) = StdExponential
@inline MeasureBase.transport_to_std(::Type{StdExponential}, d::Weibull, x) = _nan_outside(d, x, abs(x / d.θ)^d.α)
@inline MeasureBase.transport_from_std(::Type{StdExponential}, d::Weibull, z) = _nan_outside(StdExponential, z, d.θ * abs(z)^(1 / d.α))

@inline MeasureBase.preferred_stdmeasure(::Type{<:Gamma}) = StdUniform
@inline MeasureBase.transport_to_std(::Type{StdUniform}, d::Gamma, x) = _nan_outside(d, x, _gamma_cdf(d.α, _nonnan(abs(x / d.θ), 1)))
@inline MeasureBase.transport_from_std(::Type{StdUniform}, d::Gamma, p) = _nan_outside(StdUniform, p, d.θ * _gamma_quantile(d.α, _quantile_arg(p)))

@inline MeasureBase.preferred_stdmeasure(::Type{<:Beta}) = StdUniform
@inline MeasureBase.transport_to_std(::Type{StdUniform}, d::Beta, x) = _nan_outside(d, x, _beta_cdf(d.α, d.β, _nonnan(_unit_clamp(x), 1 // 2)))
@inline MeasureBase.transport_from_std(::Type{StdUniform}, d::Beta, p) = _nan_outside(StdUniform, p, _beta_quantile(d.α, d.β, _quantile_arg(p)))

@inline _unit_clamp(x) = clamp(x, zero(x), one(x))

# The Stieltjes measure functions of the gamma and beta families and their
# inverses go through the kernels, so that they are differentiable with
# respect to the shapes (e.g. in truncations). The complements use the
# complementary gamma kernels resp. the mirrored beta family to keep the
# precision of the upper tail. The kernels only see arguments inside the
# support, the values outside are exact:
@inline function MeasureBase.smf(d::Gamma, x::Real)
    z = x / d.θ
    _nan_or(x, _support_select(z > 0, z < Inf, _gamma_cdf(d.α, _gamma_arg(z)), false))
end
@inline function MeasureBase.csmf(d::Gamma, x::Real)
    z = x / d.θ
    _nan_or(x, _support_select(z > 0, z < Inf, _gamma_ccdf(d.α, _gamma_arg(z)), true))
end
@inline MeasureBase.invsmf(d::Gamma, p::Real) = _nan_outside(StdUniform, p, d.θ * _gamma_quantile(d.α, _quantile_arg(p)))
@inline MeasureBase.invcsmf(d::Gamma, p::Real) = _nan_outside(StdUniform, p, d.θ * _gamma_cquantile(d.α, _quantile_arg(p)))

@inline function MeasureBase.smf(d::Beta, x::Real)
    _nan_or(x, _support_select(x > 0, x < 1, _beta_cdf(d.α, d.β, _beta_arg(x)), false))
end
@inline function MeasureBase.csmf(d::Beta, x::Real)
    _nan_or(x, _support_select(x > 0, x < 1, _beta_cdf(d.β, d.α, 1 - _beta_arg(x)), true))
end
@inline MeasureBase.invsmf(d::Beta, p::Real) = _nan_outside(StdUniform, p, _beta_quantile(d.α, d.β, _quantile_arg(p)))
@inline MeasureBase.invcsmf(d::Beta, p::Real) = _nan_outside(StdUniform, p, 1 - _beta_quantile(d.β, d.α, _quantile_arg(p)))

@inline _gamma_arg(z) = (u = float(z); ifelse((u > 0) & (u < Inf), u, one(u)))
@inline _beta_arg(x) = (u = float(x); ifelse((u > 0) & (u < 1), u, oftype(u, 1 // 2)))
# `y` inside the support, below and above it the limits of the Stieltjes
# measure function, resp. of its complement if `compl`:
@inline _support_select(above_lo, below_hi, y, compl::Bool) =
    ifelse(above_lo, ifelse(below_hi, y, oftype(y, !compl)), oftype(y, compl))
@inline _nan_or(x, y) = ifelse(isnan(x), oftype(y, NaN), y)
