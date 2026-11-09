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
@inline MeasureBase.logdensityof_impl(m::AsMeasure{<:_Families}, x) =
    MeasureBase._checksupport(MeasureBase.insupportof(m, x), _family_logd(m.obj, x))

# Zero by value: dual numbers compare lexicographically with their
# partials, but not across the smallest normal float (subnormals count as
# zero, as on devices that flush them):
@inline _iszero_value(x) = abs(x) < floatmin(float(typeof(x)))

# `c * log(y)` for `y >= 0`, zero for `c == y == 0`. Only that endpoint is
# guarded, so that the derivative with respect to `c` survives at `c == 0`.
# Integer `y` take the float type of `c`:
@inline function _clog(c, y)
    u = y + zero(c)
    c * log(ifelse(_iszero_value(c) & _iszero_value(u), one(u), u))
end
# `c * log(1 - x)` for `x <= 1`, zero for `c == 0, x == 1`:
@inline function _clog1m(c, x)
    u = x + zero(c)
    c * log1p(-ifelse(_iszero_value(c) & (u == 1), zero(u), u))
end

# For large shape parameters the terms of the direct density formulas of
# the gamma, beta and Poisson families are large and cancel to a small
# fraction of their size. Above `_concentrated` the densities expand
# around the mode instead, with Stirling's series for `loggamma`: with
# `y = a (1 + t)`, `(a - 1) log(y / a) - (y - a) = a (log1p(t) - t) - log1p(t)`.
const _concentrated = 10
@inline _atleast_concentrated(a) = ifelse(a >= _concentrated, a, _concentrated * one(a))

# `log1p(t)` for `t = r - 1`, from `r` where `1 + t` would lose its
# relative precision:
@inline _log1p_of(t, r) = ifelse(t < -one(t) / 2, log(r), log1p(t))

# Log-density of the gamma distribution with unit scale at `y >= 0`:
@inline function _gamma_logd(α, y)
    direct = _clog(α - 1, y) - y - loggamma(α)
    a = _atleast_concentrated(α)
    inside = (y > 0) & (y < Inf)
    yc = ifelse(inside, y, a)
    t = (yc - a) / a
    lt = _log1p_of(t, yc / a)
    expanded = a * (lt - t) - lt - (log(a) + log2π) / 2 - _stirling_err(a)
    ifelse(α >= _concentrated, ifelse(inside, expanded, oftype(expanded, -Inf)), direct)
end

# Log-density of the beta distribution at `0 <= x <= 1`:
@inline function _beta_logd(α, β, x)
    direct = _clog(α - 1, x) + _clog1m(β - 1, x) - logbeta(α, β)
    a, b = _atleast_concentrated(α), _atleast_concentrated(β)
    n = a + b
    inside = (x > 0) & (x < 1)
    xc = ifelse(inside, x, one(x) / 2)
    nx, ny = n * xc, n * (1 - xc)
    u, v = (nx - a) / a, (ny - b) / b
    lu, lv = _log1p_of(u, nx / a), _log1p_of(v, ny / b)
    expanded = a * (lu - u) + b * (lv - v) - log(xc) - log1p(-xc) +
        (log(a) + log(b) - log(n) - log2π) / 2 - (_stirling_err(a) + _stirling_err(b) - _stirling_err(n))
    ifelse(min(α, β) >= _concentrated, ifelse(inside, expanded, oftype(expanded, -Inf)), direct)
end

# Log-probability of the Poisson distribution at `|k|`:
@inline function _poisson_logd(λ, k)
    kf = abs(k) + zero(float(λ))
    direct = _clog(kf, λ) - λ - loggamma(kf + 1)
    kc, λc = _atleast_concentrated(kf), _atleast_concentrated(λ)
    t = (kc - λc) / λc
    lt = _log1p_of(t, kc / λc)
    expanded = (kc - λc) - kc * lt - (log(kc) + log2π) / 2 - _stirling_err(kc)
    ifelse((kf >= _concentrated) & (λ >= _concentrated), expanded, direct)
end

# Constants stay irrational, so that they take the precision of the variates
# and parameters:
@inline function _family_logd(d::Normal, x)
    z = (x - d.μ) / d.σ
    -(z * z + log2π) / 2 - log(d.σ)
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
    ℓ = -(z * z + log2π) / 2 - log(d.σ) - lx
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
    ℓ = _gamma_logd(d.α, abs(x) / d.θ) - log(d.θ)
    ifelse(isinf(x), oftype(ℓ, -Inf), ℓ)
end
@inline _family_insupport(::Gamma, x) = x >= 0

@inline _family_logd(d::Beta, x) = _beta_logd(d.α, d.β, _unit_clamp(x))
@inline _family_insupport(::Beta, x) = (0 <= x) & (x <= 1)

@inline _family_logd(d::Poisson, x) = _poisson_logd(d.λ, x)
@inline _family_insupport(::Poisson, x) = (x >= 0) & (x == floor(x))

@inline _family_logd(d::Bernoulli, x) = ifelse(x == 1, log(d.p), log1p(-d.p))
@inline _family_insupport(::Bernoulli, x) = (x == 0) | (x == 1)


# Standard transports of the non-affine families: Cauchy, Laplace, Gamma
# and Beta pivot on the logistic measure, the log-normal and Weibull
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
@inline MeasureBase._gamma_logpdf(α, x) = _gamma_logd(α, x)
@inline MeasureBase._beta_logpdf(α, β, x) = _beta_logd(α, β, x)

# Branch-free formulas evaluate all branches, also in AD pullbacks, so
# unselected branches must stay finite with finite derivatives. Tails are
# computed from tail probabilities, selected by the same condition as the
# result, so that derivatives at the median follow one branch.

# The lower and upper tail probabilities `p` and `q` of a variate give the
# standard logistic variate `log(p) - log(q)`, which resolves both tails
# (with the probabilities floored at the smallest positive float). Back
# from it, the quantile comes from the smaller one of `logistic(±l)`:
@inline _logit_of_tails(p, q) = log(MeasureBase._positive_prob(p)) - log(MeasureBase._positive_prob(q))

# `log(p)` from `q` where `p` is close to one:
@inline function _log_tail(p, q)
    lower = p < q
    ifelse(lower, log(ifelse(lower, p, one(p))), log1p(-ifelse(lower, zero(q), q)))
end

# `|x|` with the derivative of `x` at zero (`abs` has none there in some AD
# backends):
@inline _signed_abs(x) = ifelse(x < zero(x), -x, x)

# The (log) probability arguments of the inverse SMFs in both tails:
@inline _probs_of_logp(logp) = (lp = _logquantile_arg(logp); (exp(lp), -expm1(lp)))

# Cauchy: `atan(z) / π + 1/2` cancels in the tails, the angles
# `atan(1, -z)` and `atan(1, z)` don't. `tan(π(p - 1/2))` loses small
# probabilities in `p - 1/2`, the cotangent of the smaller tail probability
# doesn't.
@inline _cauchy_cdf(d::Cauchy, x) = (z = (x - d.μ) / d.σ; atan(one(z), -z) / π)
@inline _cauchy_ccdf(d::Cauchy, x) = (z = (x - d.μ) / d.σ; atan(one(z), z) / π)
@inline function _cauchy_tail_quantile(d::Cauchy, p, q)
    lower = p < q
    r = _unit_interior(ifelse(lower, p, q))
    t = cospi(r) / sinpi(r)
    muladd(d.σ, ifelse(lower, -t, t), d.μ)
end

@inline MeasureBase.preferred_stdmeasure(::Type{<:Cauchy}) = StdLogistic
@inline function MeasureBase.transport_to_std(::Type{StdLogistic}, d::Cauchy, x)
    _nan_outside(d, x, _logit_of_tails(_cauchy_cdf(d, x), _cauchy_ccdf(d, x)))
end
@inline function MeasureBase.transport_from_std(::Type{StdLogistic}, d::Cauchy, l)
    _cauchy_tail_quantile(d, logistic(l), logistic(-l))
end

@inline MeasureBase.smf(d::Cauchy, x::Real) = _dist_smf(_cauchy_cdf, d, x, x > 0)
@inline MeasureBase.csmf(d::Cauchy, x::Real) = _dist_smf(_cauchy_ccdf, d, x, x < 0)
@inline MeasureBase.logsmf(d::Cauchy, x::Real) = _dist_logsmf((d, x) -> _log_tail(_cauchy_cdf(d, x), _cauchy_ccdf(d, x)), d, x, x > 0)
@inline MeasureBase.logcsmf(d::Cauchy, x::Real) = _dist_logsmf((d, x) -> _log_tail(_cauchy_ccdf(d, x), _cauchy_cdf(d, x)), d, x, x < 0)
@inline MeasureBase.invsmf(d::Cauchy, p::Real) = _nan_outside(StdUniform, p, _cauchy_tail_quantile(d, p, 1 - p))
@inline MeasureBase.invcsmf(d::Cauchy, q::Real) = _nan_outside(StdUniform, q, _cauchy_tail_quantile(d, 1 - q, q))
@inline MeasureBase.invlogsmf(d::Cauchy, logp::Real) = _nan_above_zero(logp, _cauchy_tail_quantile(d, _probs_of_logp(logp)...))
@inline function MeasureBase.invlogcsmf(d::Cauchy, logq::Real)
    q, p = _probs_of_logp(logq)
    _nan_above_zero(logq, _cauchy_tail_quantile(d, p, q))
end

# Laplace: `l = ±log(2 exp(|z|) - 1)` by the sign of `z`, and back
# `|z| = log((exp(|l|) + 1) / 2)`, each in two forms that stay accurate
# near zero resp. don't overflow:
@inline MeasureBase.preferred_stdmeasure(::Type{<:Laplace}) = StdLogistic
@inline function MeasureBase.transport_to_std(::Type{StdLogistic}, d::Laplace, x)
    z = (x - d.μ) / d.θ
    a = _signed_abs(z)
    g = ifelse(a < 1, log1p(2 * expm1(min(a, one(a)))), a + logtwo + log1p(-exp(-a) / 2))
    _nan_outside(d, x, ifelse(z < 0, -g, g))
end
@inline function MeasureBase.transport_from_std(::Type{StdLogistic}, d::Laplace, l)
    b = _signed_abs(l)
    a = ifelse(b < 1, log1p(expm1(min(b, one(b))) / 2), log1pexp(b) - logtwo)
    muladd(d.θ, ifelse(l < 0, -a, a), d.μ)
end

@inline MeasureBase.preferred_stdmeasure(::Type{<:LogNormal}) = StdNormal
# The lower end of the support maps to a finite standard variate:
@inline MeasureBase.transport_to_std(::Type{StdNormal}, d::LogNormal, x) =
    _nan_outside(d, x, (log(MeasureBase._positive_prob(abs(x))) - d.μ) / d.σ)
@inline MeasureBase.transport_from_std(::Type{StdNormal}, d::LogNormal, z) = exp(muladd(d.σ, z, d.μ))

@inline MeasureBase.preferred_stdmeasure(::Type{<:Weibull}) = StdExponential
@inline MeasureBase.transport_to_std(::Type{StdExponential}, d::Weibull, x) = _nan_outside(d, x, abs(x / d.θ)^d.α)
@inline MeasureBase.transport_from_std(::Type{StdExponential}, d::Weibull, z) = _nan_outside(StdExponential, z, d.θ * abs(z)^(1 / d.α))

# Gamma and Beta: the tail probabilities from the incomplete gamma resp.
# the (mirrored) incomplete beta functions, quantiles from the smaller one:
@inline function _gamma_tail_quantile(α, p, q)
    lower = p < q
    r = _quantile_arg(ifelse(lower, p, q))
    h = one(r) / 2
    ifelse(lower, _gamma_quantile(α, ifelse(lower, r, h)), _gamma_cquantile(α, ifelse(lower, h, r)))
end
# Beta quantiles from both tail probabilities, the variate `x` from the
# lower and its complement `1 - x` from the mirrored one, the smaller one
# of them is accurate:
@inline function _beta_quantiles(α, β, p, q)
    x, y = _beta_quantile(α, β, _quantile_arg(p)), _beta_quantile(β, α, _quantile_arg(q))
    small_x = x <= y
    (ifelse(small_x, x, 1 - y), ifelse(small_x, 1 - x, y))
end
@inline _beta_tail_quantile(α, β, p, q) = first(_beta_quantiles(α, β, p, q))

# Both tail probabilities of Beta(α, β) at `z` with the complement
# `c = 1 - z`, each given to relative precision. The incomplete beta
# function of the smaller one of `z` and `c` is accurate, the other tail
# probability comes from the mirrored function where the rounding of its
# argument costs less than forming the complement (the density times the
# argument below one):
@inline function _beta_tails(α, β, z, c)
    p_direct, q_direct = _beta_cdf(α, β, z), _beta_cdf(β, α, c)
    ℓ = _beta_logd(α, β, z)
    p = ifelse((z > c) & !(ℓ + log(z) < 0), 1 - q_direct, p_direct)
    q = ifelse((z <= c) & !(ℓ + log(c) < 0), 1 - p_direct, q_direct)
    (p, q)
end

@inline MeasureBase.preferred_stdmeasure(::Type{<:Gamma}) = StdLogistic
@inline function MeasureBase.transport_to_std(::Type{StdLogistic}, d::Gamma, x)
    u = _gamma_trafo_arg(x / d.θ)
    _nan_outside(d, x, _logit_of_tails(_gamma_cdf(d.α, u), _gamma_ccdf(d.α, u)))
end
@inline function MeasureBase.transport_from_std(::Type{StdLogistic}, d::Gamma, l)
    _nan_or(l, d.θ * _gamma_tail_quantile(d.α, logistic(l), logistic(-l)))
end

@inline MeasureBase.preferred_stdmeasure(::Type{<:Beta}) = StdLogistic
@inline function MeasureBase.transport_to_std(::Type{StdLogistic}, d::Beta, x)
    u = _nonnan(_unit_clamp(x), 1 // 2)
    _nan_outside(d, x, _logit_of_tails(_beta_tails(d.α, d.β, u, 1 - u)...))
end
@inline function MeasureBase.transport_from_std(::Type{StdLogistic}, d::Beta, l)
    _nan_or(l, _beta_tail_quantile(d.α, d.β, logistic(l), logistic(-l)))
end

@inline _unit_clamp(x) = clamp(x, zero(x), one(x))
# `|z|` inside the support and on its lower end, one for infinite and NaN `z`:
@inline _gamma_trafo_arg(z) = ifelse(isfinite(z), abs(z), one(z))

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
@inline function MeasureBase.logsmf(d::Gamma, x::Real)
    z = x / d.θ
    u = _gamma_arg(z)
    _nan_or(x, _log_support_select(z > 0, z < Inf, _log_tail(_gamma_cdf(d.α, u), _gamma_ccdf(d.α, u)), false))
end
@inline function MeasureBase.logcsmf(d::Gamma, x::Real)
    z = x / d.θ
    u = _gamma_arg(z)
    _nan_or(x, _log_support_select(z > 0, z < Inf, _log_tail(_gamma_ccdf(d.α, u), _gamma_cdf(d.α, u)), true))
end
@inline MeasureBase.invsmf(d::Gamma, p::Real) = _nan_outside(StdUniform, p, d.θ * _gamma_quantile(d.α, _quantile_arg(p)))
@inline MeasureBase.invcsmf(d::Gamma, p::Real) = _nan_outside(StdUniform, p, d.θ * _gamma_cquantile(d.α, _quantile_arg(p)))
@inline MeasureBase.invlogsmf(d::Gamma, logp::Real) = _nan_above_zero(logp, d.θ * _gamma_tail_quantile(d.α, _probs_of_logp(logp)...))
@inline function MeasureBase.invlogcsmf(d::Gamma, logq::Real)
    q, p = _probs_of_logp(logq)
    _nan_above_zero(logq, d.θ * _gamma_tail_quantile(d.α, p, q))
end

@inline function MeasureBase.smf(d::Beta, x::Real)
    _nan_or(x, _support_select(x > 0, x < 1, _beta_cdf(d.α, d.β, _beta_arg(x)), false))
end
@inline function MeasureBase.csmf(d::Beta, x::Real)
    u = _beta_arg(x)
    _nan_or(x, _support_select(x > 0, x < 1, last(_beta_tails(d.α, d.β, u, 1 - u)), true))
end
@inline function MeasureBase.logsmf(d::Beta, x::Real)
    u = _beta_arg(x)
    _nan_or(x, _log_support_select(x > 0, x < 1, _log_tail(_beta_tails(d.α, d.β, u, 1 - u)...), false))
end
@inline function MeasureBase.logcsmf(d::Beta, x::Real)
    u = _beta_arg(x)
    p, q = _beta_tails(d.α, d.β, u, 1 - u)
    _nan_or(x, _log_support_select(x > 0, x < 1, _log_tail(q, p), true))
end
@inline MeasureBase.invsmf(d::Beta, p::Real) = _nan_outside(StdUniform, p, _beta_tail_quantile(d.α, d.β, p, 1 - p))
@inline MeasureBase.invcsmf(d::Beta, q::Real) = _nan_outside(StdUniform, q, _beta_tail_quantile(d.α, d.β, 1 - q, q))
@inline MeasureBase.invlogsmf(d::Beta, logp::Real) = _nan_above_zero(logp, _beta_tail_quantile(d.α, d.β, _probs_of_logp(logp)...))
@inline function MeasureBase.invlogcsmf(d::Beta, logq::Real)
    q, p = _probs_of_logp(logq)
    _nan_above_zero(logq, _beta_tail_quantile(d.α, d.β, p, q))
end

@inline _gamma_arg(z) = (u = float(z); ifelse((u > 0) & (u < Inf), u, one(u)))
@inline _beta_arg(x) = (u = float(x); ifelse((u > 0) & (u < 1), u, oftype(u, 1 // 2)))
# `y` inside the support, below and above it the limits of the Stieltjes
# measure function, resp. of its complement if `compl`, and of their logs:
@inline _support_select(above_lo, below_hi, y, compl::Bool) =
    ifelse(above_lo, ifelse(below_hi, y, oftype(y, !compl)), oftype(y, compl))
@inline _log_support_select(above_lo, below_hi, y, compl::Bool) =
    ifelse(above_lo, ifelse(below_hi, y, oftype(y, compl ? -Inf : 0)), oftype(y, compl ? 0 : -Inf))
@inline _nan_or(x, y) = ifelse(isnan(x), oftype(y, NaN), y)


# Families of gamma and beta variates: chi-squared and Erlang distributions
# are gamma distributions, inverse gamma variates are reciprocals of gamma
# variates, beta prime and F variates are ratios `(b / a) y / (1 - y)` of
# beta variates `y`. Their Stieltjes measure functions, quantiles and
# standard transports go through the gamma and beta kernels, so that they
# are differentiable with respect to the shapes in all AD backends:

const _GammaReparam = Union{Distributions.Chisq,Distributions.Erlang}

@inline _as_gamma(d::Distributions.Chisq) = Gamma(Distributions.dof(d) / 2, 2; check_args = false)
@inline _as_gamma(d::Distributions.Erlang) = Gamma(Distributions.shape(d), Distributions.scale(d); check_args = false)

for f in (:smf, :csmf, :logsmf, :logcsmf, :invsmf, :invcsmf, :invlogsmf, :invlogcsmf)
    @eval @inline MeasureBase.$f(d::_GammaReparam, x::Real) = MeasureBase.$f(_as_gamma(d), x)
end
@inline MeasureBase.transport_to_std(::Type{StdLogistic}, d::_GammaReparam, x) = transport_to_std(StdLogistic, _as_gamma(d), x)
@inline MeasureBase.transport_from_std(::Type{StdLogistic}, d::_GammaReparam, l) = transport_from_std(StdLogistic, _as_gamma(d), l)

# Quantiles that overflow stay at the largest finite value:
@inline _capped(x) = min(x, floatmax(typeof(x)))

# Inverse gamma variates `x = θ / z` of standard gamma variates `z`, the
# lower tail of `x` is the upper one of `z`:
@inline _std_gamma(d::Distributions.InverseGamma) = Gamma(Distributions.shape(d), 1; check_args = false)
@inline _invgamma_arg(d::Distributions.InverseGamma, x) = Distributions.scale(d) / ifelse(x > 0, x, one(x))

@inline function MeasureBase.smf(d::Distributions.InverseGamma, x::Real)
    F = MeasureBase.csmf(_std_gamma(d), _invgamma_arg(d, x))
    _nan_or(x, ifelse(x > 0, F, zero(F)))
end
@inline function MeasureBase.csmf(d::Distributions.InverseGamma, x::Real)
    G = MeasureBase.smf(_std_gamma(d), _invgamma_arg(d, x))
    _nan_or(x, ifelse(x > 0, G, one(G)))
end
@inline function MeasureBase.logsmf(d::Distributions.InverseGamma, x::Real)
    ℓ = MeasureBase.logcsmf(_std_gamma(d), _invgamma_arg(d, x))
    _nan_or(x, ifelse(x > 0, ℓ, oftype(ℓ, -Inf)))
end
@inline function MeasureBase.logcsmf(d::Distributions.InverseGamma, x::Real)
    ℓ = MeasureBase.logsmf(_std_gamma(d), _invgamma_arg(d, x))
    _nan_or(x, ifelse(x > 0, ℓ, zero(ℓ)))
end
for (f, f_gamma) in ((:invsmf, :invcsmf), (:invcsmf, :invsmf), (:invlogsmf, :invlogcsmf), (:invlogcsmf, :invlogsmf))
    @eval @inline MeasureBase.$f(d::Distributions.InverseGamma, p::Real) =
        _capped(Distributions.scale(d) / MeasureBase.$f_gamma(_std_gamma(d), p))
end

@inline function MeasureBase.transport_to_std(::Type{StdLogistic}, d::Distributions.InverseGamma, x)
    _nan_outside(d, x, _logit_of_tails(MeasureBase.smf(d, x), MeasureBase.csmf(d, x)))
end
@inline function MeasureBase.transport_from_std(::Type{StdLogistic}, d::Distributions.InverseGamma, l)
    z = _gamma_tail_quantile(Distributions.shape(d), logistic(-l), logistic(l))
    _nan_or(l, _capped(Distributions.scale(d) / z))
end

# The beta shapes `(α, β)` and the weights `(a, b)` of the ratio:
const _BetaRatio = Union{Distributions.BetaPrime,Distributions.FDist}
@inline _beta_ratio_params(d::Distributions.BetaPrime) = (Distributions.params(d)..., 1, 1)
@inline _beta_ratio_params(d::Distributions.FDist) = ((ν1, ν2) = Distributions.params(d); (ν1 / 2, ν2 / 2, ν1, ν2))

# Both tail probabilities at `x`, from the beta variate `y = r / (1 + r)`
# and its complement `1 / (1 + r)` for `r = a x / b`, each formed without
# cancellation. Exact limits outside the support:
@inline function _beta_ratio_tails(d::_BetaRatio, x)
    α, β, a, b = _beta_ratio_params(d)
    inside = (x > 0) & (x < Inf)
    r = a * ifelse(inside, x, one(x)) / b
    big = r > 1
    t = ifelse(big, inv(ifelse(big, r, one(r))), r)
    y, c = ifelse(big, one(t), t) / (1 + t), ifelse(big, t, one(t)) / (1 + t)
    p, q = _beta_tails(α, β, y, c)
    (_support_select(x > 0, x < Inf, p, false), _support_select(x > 0, x < Inf, q, true))
end

@inline function _beta_ratio_quantile(d::_BetaRatio, p, q)
    α, β, a, b = _beta_ratio_params(d)
    y, c = _beta_quantiles(α, β, p, q)
    _capped(b * y / (a * c))
end

@inline MeasureBase.smf(d::_BetaRatio, x::Real) = _nan_or(x, first(_beta_ratio_tails(d, x)))
@inline MeasureBase.csmf(d::_BetaRatio, x::Real) = _nan_or(x, last(_beta_ratio_tails(d, x)))
@inline function MeasureBase.logsmf(d::_BetaRatio, x::Real)
    p, q = _beta_ratio_tails(d, x)
    _nan_or(x, _log_support_select(x > 0, x < Inf, _log_tail(p, q), false))
end
@inline function MeasureBase.logcsmf(d::_BetaRatio, x::Real)
    p, q = _beta_ratio_tails(d, x)
    _nan_or(x, _log_support_select(x > 0, x < Inf, _log_tail(q, p), true))
end
@inline MeasureBase.invsmf(d::_BetaRatio, p::Real) = _nan_outside(StdUniform, p, _beta_ratio_quantile(d, p, 1 - p))
@inline MeasureBase.invcsmf(d::_BetaRatio, q::Real) = _nan_outside(StdUniform, q, _beta_ratio_quantile(d, 1 - q, q))
@inline MeasureBase.invlogsmf(d::_BetaRatio, logp::Real) = _nan_above_zero(logp, _beta_ratio_quantile(d, _probs_of_logp(logp)...))
@inline function MeasureBase.invlogcsmf(d::_BetaRatio, logq::Real)
    q, p = _probs_of_logp(logq)
    _nan_above_zero(logq, _beta_ratio_quantile(d, p, q))
end

@inline MeasureBase.transport_to_std(::Type{StdLogistic}, d::_BetaRatio, x) =
    _nan_outside(d, x, _logit_of_tails(_beta_ratio_tails(d, x)...))
@inline MeasureBase.transport_from_std(::Type{StdLogistic}, d::_BetaRatio, l) =
    _nan_or(l, _beta_ratio_quantile(d, logistic(l), logistic(-l)))
