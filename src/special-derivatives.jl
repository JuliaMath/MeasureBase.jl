# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using SpecialFunctions: digamma, loggamma, logbeta
using LogExpFunctions: log1pmx

# Derivatives of the regularized incomplete gamma and beta functions with
# respect to their shape parameters. The power series and continued
# fractions of the functions are differentiated term by term, following
#
# * R. J. Moore, "Algorithm AS 187: Derivatives of the Incomplete Gamma
#   Integral", J. R. Stat. Soc. C 31, 330 (1982)
# * R. J. Boik and J. F. Robison-Cox, "Derivatives of the Incomplete Beta
#   Function", J. Stat. Softw. 3, 1 (1998)
#
# The series and continued fractions are summed until they converge. Near
# the median they need a number of terms that grows with the square root of
# the shape parameters, and the derivatives are NaN where the series don't
# converge within `_special_maxiter` terms (shape parameters above about
# 10^12). Like the incomplete gamma and beta functions of SpecialFunctions,
# these run on the host, not on GPUs or under Reactant.

const _special_maxiter = 10_000_000

# Near the median, the series and continued fractions need up to about
# 9 sqrt(shape) terms:
_special_budget(shape::Real) = trunc(Int, min(100_000 + 20 * sqrt(shape), _special_maxiter))

# log(r) for r = 1 + d, from r where 1 + d has lost relative precision:
_log1p_or_log(d::Real, r::Real) = d < -one(d) / 2 ? log(r) : log1p(d)

# Asymptotic expansions of log(z) - ψ(z) and of the Stirling error
# log Γ(z) - ((z - 1/2) log(z) - z + log(2π) / 2), for z ≥ 10:
function _log_minus_digamma(z::T) where {T<:Real}
    u = inv(z)^2
    inv(2z) + u * evalpoly(u, map(T, (1//12, -1//120, 1//252, -1//240, 1//132, -691//32760, 1//12)))
end
function _stirling_err(z::T) where {T<:Real}
    u = inv(z)^2
    evalpoly(u, map(T, (1//12, -1//360, 1//1260, -1//1680, 1//1188, -691//360360, 1//156))) / z
end

# The derivative ∂P(α, x)/∂α of the regularized lower incomplete gamma
# function P(α, x):
_gamma_cdf_dshape(α::Real, x::Real) = _gamma_cdf_dshape_impl(promote(float(α), float(x))...)

_gamma_cdf_dshape_impl(α::T, x::T) where {T<:Union{Float16,Float32}} = T(_gamma_cdf_dshape_impl(Float64(α), Float64(x)))

function _gamma_cdf_dshape_impl(α::T, x::T) where {T<:Real}
    valid_α = (α > 0) & (α < Inf)
    valid_α & (x > 0) & (x < Inf) || return ifelse(valid_α & ((x == 0) | (x == Inf)), zero(T), convert(T, NaN))
    (x > 1) & (x > α) ? -_gamma_ccdf_dshape_cf(α, x) : _gamma_cdf_dshape_series(α, x)
end

# log(x^α e^(-x) / Γ(α + k)) and log(x) - ψ(α + k) for k = 0, 1, for large
# α without the cancellation of the terms of order α log(α):
function _gamma_log_terms(α::T, x::T, k::Integer) where {T<:Real}
    if α >= 10
        t = (x - α) / α
        l = _log1p_or_log(t, x / α)
        lmt = t < -one(t) / 2 ? l - t : log1pmx(t)
        logK = α * lmt + log(α / (2 * T(π))) / 2 - _stirling_err(α) - k * log(α)
        logK, l + _log_minus_digamma(α) - k / α
    else
        α * log(x) - x - loggamma(α + k), log(x) - digamma(α + k)
    end
end

# P(α, x) = x^α e^(-x) / Γ(α + 1) * ∑ₙ cₙ with c₀ = 1, cₙ = cₙ₋₁ x / (α + n):
function _gamma_cdf_dshape_series(α::T, x::T) where {T<:Real}
    tol = eps(T)
    logK, L = _gamma_log_terms(α, x, 1)
    c, dc = one(T), zero(T)
    s, ds = c, dc
    for n in 1:_special_budget(α)
        r = x / (α + n)
        dc = r * (dc - c / (α + n))
        c *= r
        s += c
        ds += dc
        if (c <= tol * s) & (abs(dc) <= tol * (abs(ds) + s * abs(L)))
            return exp(logK) * (s * L + ds)
        end
    end
    return convert(T, NaN)
end

# Q(α, x) = 1 - P(α, x) = x^α e^(-x) / Γ(α) * F with the continued fraction
# F = 1 / (b₀ + a₁ / (b₁ + a₂ / (b₂ + ...))), aₖ = k (α - k), bₖ = x + 2k + 1 - α:
function _gamma_ccdf_dshape_cf(α::T, x::T) where {T<:Real}
    tol = eps(T)
    logK, L = _gamma_log_terms(α, x, 0)
    # Convergents Aₙ / Bₙ and their derivatives by the Wallis recurrences,
    # normalized to Bₙ = 1:
    A₀, B₀, dA₀, dB₀ = one(T), zero(T), zero(T), zero(T)
    A₁, B₁, dA₁, dB₁ = zero(T), one(T), zero(T), zero(T)
    F, dF = zero(T), zero(T)
    for k in 0:_special_budget(α)
        a, da = ifelse(k == 0, one(T), k * (α - k)), convert(T, k)
        b = x + (2k + 1) - α
        A₂, B₂ = b * A₁ + a * A₀, b * B₁ + a * B₀
        dA₂ = b * dA₁ - A₁ + a * dA₀ + da * A₀
        dB₂ = b * dB₁ - B₁ + a * dB₀ + da * B₀
        r = inv(B₂)
        A₀, B₀, dA₀, dB₀ = A₁ * r, B₁ * r, dA₁ * r, dB₁ * r
        A₁, B₁, dA₁, dB₁ = A₂ * r, one(T), dA₂ * r, dB₂ * r
        F_prev, dF_prev = F, dF
        F, dF = A₁, dA₁ - A₁ * dB₁
        if (abs(F - F_prev) <= tol * abs(F)) & (abs(dF - dF_prev) <= tol * (abs(dF) + abs(F * L)))
            return exp(logK) * (F * L + dF)
        end
    end
    return convert(T, NaN)
end


# The derivatives (∂I_x(α, β)/∂α, ∂I_x(α, β)/∂β) of the regularized
# incomplete beta function I_x(α, β):
_beta_cdf_dshape(α::Real, β::Real, x::Real) = _beta_cdf_dshape_impl(promote(float(α), float(β), float(x))...)

function _beta_cdf_dshape_impl(α::T, β::T, x::T) where {T<:Union{Float16,Float32}}
    map(T, _beta_cdf_dshape_impl(Float64(α), Float64(β), Float64(x)))
end

function _beta_cdf_dshape_impl(α::T, β::T, x::T) where {T<:Real}
    valid_shapes = (α > 0) & (β > 0) & (α < Inf) & (β < Inf)
    if !(valid_shapes & (x > 0) & (x < 1))
        d = ifelse(valid_shapes & ((x == 0) | (x == 1)), zero(T), convert(T, NaN))
        return (d, d)
    end
    y = 1 - x
    logS = _beta_log_prefactor(α, β, x)
    # The continued fraction converges fast below (α + 1) / (α + β + 2), use
    # I_x(α, β) = 1 - I_{1-x}(β, α) above:
    if x * (α + β + 2) < α + 1
        _beta_cdf_dshape_cf(α, β, x, log(x), log1p(-x), logS)
    else
        dβ, dα = _beta_cdf_dshape_cf(β, α, y, log(y), log(x), logS)
        (-dα, -dβ)
    end
end

# log(x^α (1 - x)^β / B(α, β)), for large α and β without the cancellation
# of the terms of order α and β. With p = α / (α + β), the first order
# terms in x - p cancel, so the rounding of p doesn't matter:
function _beta_log_prefactor(α::T, β::T, x::T) where {T<:Real}
    if min(α, β) >= 10
        s = α + β
        p, q = α / s, β / s
        δ = x - p
        α * _log1p_or_log(δ / p, x / p) + β * _log1p_or_log(-δ / q, (1 - x) / q) +
            log(α * q / (2 * T(π))) / 2 - _stirling_err(α) - _stirling_err(β) + _stirling_err(s)
    else
        α * log(x) + β * log1p(-x) - logbeta(α, β)
    end
end

# I_x(α, β) = x^α y^β / (α B(α, β)) * F with y = 1 - x and the continued
# fraction F = 1 / (1 + d₁ / (1 + d₂ / (1 + ...))),
# d₂ₘ = m (β - m) x / ((α + 2m - 1) (α + 2m)),
# d₂ₘ₊₁ = -(α + m) (α + β + m) x / ((α + 2m) (α + 2m + 1)), where
# logS = log(x^α y^β / B(α, β)):
function _beta_cdf_dshape_cf(α::T, β::T, x::T, logx::T, logy::T, logS::T) where {T<:Real}
    tol = eps(T)
    logK = logS - log(α)
    Lα = logx + _digamma_diff(α + 1, β - 1)
    Lβ = logy + _digamma_diff(β, α)
    # Convergents Aₙ / Bₙ and their derivatives by the Wallis recurrences,
    # normalized to Bₙ = 1:
    A₀, B₀, A₁, B₁ = zero(T), one(T), one(T), one(T)
    dαA₀, dαB₀, dαA₁, dαB₁ = zero(T), zero(T), zero(T), zero(T)
    dβA₀, dβB₀, dβA₁, dβB₁ = zero(T), zero(T), zero(T), zero(T)
    F, dαF, dβF = one(T), zero(T), zero(T)
    for n in 1:_special_budget(α + β)
        m = n >> 1
        d, dαd, dβd = if iseven(n)
            d = m * (β - m) * x / ((α + 2m - 1) * (α + 2m))
            (d, -d * (inv(α + 2m - 1) + inv(α + 2m)), m * x / ((α + 2m - 1) * (α + 2m)))
        else
            d = -(α + m) * (α + β + m) * x / ((α + 2m) * (α + 2m + 1))
            (d, d * (m / ((α + m) * (α + 2m)) + inv(α + β + m) - inv(α + 2m + 1)), d / (α + β + m))
        end
        A₂, B₂ = A₁ + d * A₀, B₁ + d * B₀
        dαA₂, dαB₂ = dαA₁ + d * dαA₀ + dαd * A₀, dαB₁ + d * dαB₀ + dαd * B₀
        dβA₂, dβB₂ = dβA₁ + d * dβA₀ + dβd * A₀, dβB₁ + d * dβB₀ + dβd * B₀
        r = inv(B₂)
        A₀, B₀, dαA₀, dαB₀, dβA₀, dβB₀ = A₁ * r, B₁ * r, dαA₁ * r, dαB₁ * r, dβA₁ * r, dβB₁ * r
        A₁, B₁, dαA₁, dαB₁, dβA₁, dβB₁ = A₂ * r, one(T), dαA₂ * r, dαB₂ * r, dβA₂ * r, dβB₂ * r
        F_prev, dαF_prev, dβF_prev = F, dαF, dβF
        F, dαF, dβF = A₁, dαA₁ - A₁ * dαB₁, dβA₁ - A₁ * dβB₁
        if (abs(F - F_prev) <= tol * F) & (abs(dαF - dαF_prev) <= tol * (abs(dαF) + F * abs(Lα))) &
            (abs(dβF - dβF_prev) <= tol * (abs(dβF) + F * abs(Lβ)))
            K = exp(logK)
            return (K * (F * Lα + dαF), K * (F * Lβ + dβF))
        end
    end
    return (convert(T, NaN), convert(T, NaN))
end

# ψ(z + a) - ψ(z) for z > 0 and z + a > 0, also for a ≪ z. Shifts z to
# z ≥ 10 by ψ(z + 1) = ψ(z) + 1/z and uses the asymptotic expansion
# ψ(z) ~ log(z) - 1/(2z) - ∑ₖ B₂ₖ / (2k z²ᵏ) there:
function _digamma_diff(z::T, a::T) where {T<:Real}
    s = zero(T)
    while z < 10
        s += a / (z * (z + a))
        z += 1
    end
    t = log1p(a / z)
    u = inv(z)^2
    # (z + a)^(-2k) - z^(-2k) = z^(-2k) expm1(-2k t):
    c = (-1 // 12, 1 // 120, -1 // 252, 1 // 240, -1 // 132, 691 // 32760, -1 // 12)
    series = sum(ntuple(k -> c[k] * u^k * expm1(-2k * t), Val(length(c))))
    s + t + a / (2 * z * (z + a)) + series
end
