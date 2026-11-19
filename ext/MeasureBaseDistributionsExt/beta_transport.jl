# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# DLMF 8.17.22: Gauss's continued fraction, evaluated by modified Lentz.
# Scalar state keeps the same evaluator usable in CPU and device kernels.
@inline function _beta_fraction(a::T, b::T, x::T, y::T) where {T}
    tiny = floatmin(T) / eps(T)
    c = one(T)
    d = y + ((one(T) - b) / (a + one(T))) * x
    d = inv(ifelse(abs(d) < tiny, copysign(tiny, d), d))
    h = d
    # Near the mean of a concentrated law, many more convergents are
    # needed. Elsewhere the residual test terminates the fraction early.
    for m in 1:100_000
        k = T(m)
        r = a + 2k
        even = (k / r) * ((b - k) / (r - one(T))) * x
        odd = -((a + k) / r) * ((a + b + k) / (r + one(T))) * x
        for coefficient in (even, odd)
            d = muladd(coefficient, d, one(T))
            c = muladd(coefficient, inv(c), one(T))
            d = inv(ifelse(abs(d) < tiny, copysign(tiny, d), d))
            c = ifelse(abs(c) < tiny, copysign(tiny, c), c)
            h *= d * c
        end
        abs(d * c - one(T)) <= 4eps(T) && return h
    end
    error("Beta transport continued fraction did not converge")
end

@inline function _beta_logprefactor(a, b, x, y, logbeta_ab)
    # The smaller coordinate survives when its complement rounds to one.
    if min(a, b) >= 10
        return if x <= y
            MeasureBase._beta_log_prefactor(a, b, x)
        else
            MeasureBase._beta_log_prefactor(b, a, y)
        end
    end
    logx, logy = x <= y ? (log(x), log1p(-x)) : (log1p(-y), log(y))
    a * logx + b * logy - logbeta_ab
end

# Near the mean of a highly concentrated law the continued fraction also
# loses precision in its first denominator. Use the standardized cumulants
# through fifth order in the Edgeworth CDF there. The first omitted terms
# are of order min(a, b)^(-2) for bounded standardized coordinates.
@inline function _beta_central_logcdf(a, b, x, y)
    s = a + b
    # Recover the rounding error in s for the centered coordinate s*x - a.
    e = min(a, b) - (s - max(a, b))
    centered = x <= y ? fma(x, s, -a) + x * e : -(fma(y, s, -b) + y * e)
    scale = sqrt(a * (b / (s + 1)))
    z = centered / scale
    abs(z) <= 4 || return nothing
    h = (b - a) / scale
    k3 = 2h / (s + 2)
    k4 = 3 * (h * k3 - 2) / (s + 3)
    k5 = 4 * (h * k4 - 6k3) / (s + 4)
    H2 = z * z - 1
    H3 = z * H2 - 2z
    H4 = z * H3 - 3H2
    H5 = z * H4 - 4H3
    H6 = z * H5 - 5H4
    H8 = z * (z * H6 - 6H5) - 7H6
    correction = k3 * H2 / 6 + k4 * H3 / 24 + k3^2 * H5 / 72 +
        k5 * H4 / 120 + k3 * k4 * H6 / 144 + k3^3 * H8 / 1296
    log(MeasureBase.Φ(z) - exp(-z^2 / 2) * invsqrt2π * correction)
end

@inline function _beta_logcdf(a, b, x, y, prefactor)
    x == 0 && return oftype(x, -Inf)
    y == 0 && return zero(y)
    if min(a, b) >= 1e8
        central = _beta_central_logcdf(a, b, x, y)
        isnothing(central) || return central
    end
    if x <= (a + one(a)) / (a + b + 2one(a))
        return min(zero(x), prefactor + log(_beta_fraction(a, b, x, y) / a))
    end
    upper = min(zero(x), prefactor + log(_beta_fraction(b, a, y, x) / b))
    log(-expm1(upper))
end

@inline function _portable_beta_cdf(a::T, b::T, x::T) where {T}
    a == b && x == T(0.5) && return x
    if T === Float32 && min(a, b) >= 1f5
        return T(_portable_beta_cdf(Float64(a), Float64(b), Float64(x)))
    end
    a == 1 && return -expm1(b * log1p(-x))
    b == 1 && return exp(a * log(x))
    y = one(T) - x
    prefactor = _beta_logprefactor(a, b, x, y, logbeta(a, b))
    exp(_beta_logcdf(a, b, x, y, prefactor))
end

@inline function _portable_beta_quantiles(a::T, b::T, p::T) where {T}
    p == 0 && return (zero(T), one(T))
    p == 1 && return (one(T), zero(T))
    a == b && p == T(0.5) && return (p, p)
    if a == 1
        logy = log1p(-p) / b
        return (-expm1(logy), exp(logy))
    elseif b == 1
        logx = log(p) / a
        return (exp(logx), -expm1(logx))
    end
    # A few ulps of CDF error across a sparse law's flat interior can move
    # its quantile by many grid points. For concentrated laws use the same
    # working precision as the forward CDF.
    if T === Float32 && ((a < 1 && b < 1) || min(a, b) >= 1f5)
        x, y = _portable_beta_quantiles(Float64(a), Float64(b), Float64(p))
        return (T(x), T(y))
    end
    reflected = p > one(T) / 2
    a, b, p = reflected ? (b, a, one(T) - p) : (a, b, p)
    logp = log(p)
    logbeta_ab = logbeta(a, b)
    bound = -log(nextfloat(zero(T)))
    # At half the smallest subnormal the power-series correction is negligible.
    cutoff = -bound - log(T(2))
    logp < a * cutoff - log(a) - logbeta_ab &&
        return reflected ? (one(T), zero(T)) : (zero(T), one(T))
    log1p(-p) < b * cutoff - log(b) - logbeta_ab &&
        return reflected ? (zero(T), one(T)) : (one(T), zero(T))
    lo, hi = -bound, bound
    mean = a / (a + b)
    initial = if min(a, b) > 1
        mean + sqrt(mean * (one(T) - mean) / (a + b + one(T))) * MeasureBase.Φinv(p)
    else
        exp((logp + log(a) + logbeta_ab) / a)
    end
    t = 0 < initial < 1 ? log(initial) - log1p(-initial) : log(a) - log(b)
    t = clamp(t, lo, hi)
    for _ in 1:128
        x, y = exp(-log1pexp(-t)), exp(-log1pexp(t))
        prefactor = _beta_logprefactor(a, b, x, y, logbeta_ab)
        value = _beta_logcdf(a, b, x, y, prefactor)
        residual = value - logp
        slope = exp(prefactor - value)
        step = residual / slope
        if abs(step) <= 8eps(T) * max(one(T), abs(t))
            return reflected ? (y, x) : (x, y)
        end
        lo, hi = residual < 0 ? (t, hi) : (lo, t)
        # Halley's correction uses d(log(S/I))/dt = a*y - b*x - S/I.
        correction = one(T) - step * (a * y - b * x - slope) / 2
        candidate = t - (abs(correction - one(T)) < one(T) / 2 ? step / correction : step)
        next = lo < candidate < hi ? candidate : (lo + hi) / 2
        if next == t
            return reflected ? (y, x) : (x, y)
        end
        t = next
    end
    error("Beta transport inverse did not converge")
end

@inline function MeasureBase._beta_cdf_impl(
    ::Type{Nothing},
    a::Real,
    b::Real,
    x::Union{Float32,Float64},
)
    aa, bb, xx = promote(float(a), float(b), x)
    _portable_beta_cdf(aa, bb, xx)
end

@inline function MeasureBase._beta_quantile_pair_impl(
    ::Type{Nothing},
    a::Real,
    b::Real,
    p::Union{Float32,Float64},
)
    aa, bb, pp = promote(float(a), float(b), p)
    _portable_beta_quantiles(aa, bb, pp)
end
