# Moments of the normalized laws of measures, where they are known exactly.

"""
    mean(μ::AbstractMeasure)

The mean of the normalized law of `μ`.

Moments are known exactly for wrapped distributions, standard measures,
point masses, products and powers (in the layout of their variates), and
truncated normal measures. Weighted and normalized measures inherit known
moments, as do affine pushforwards when AffineMaps is loaded. Other measures throw an
`ArgumentError`, MeasureBase does not estimate moments by sampling.

See also [`var`](@ref Statistics.var(::AbstractMeasure)),
[`std`](@ref Statistics.std(::AbstractMeasure)) and
[`cov`](@ref Statistics.cov(::AbstractMeasure)).
"""
Statistics.mean(μ::AbstractMeasure) = _no_moment(μ, :mean)

"""
    var(μ::AbstractMeasure)

The variance of the normalized law of `μ`, elementwise for array variates
and componentwise for tuple variates, see
[`mean`](@ref Statistics.mean(::AbstractMeasure)).
"""
Statistics.var(μ::AbstractMeasure) = _no_moment(μ, :variance)

"""
    std(μ::AbstractMeasure)

The standard deviation of the normalized law of `μ`, the square root of
[`var(μ)`](@ref Statistics.var(::AbstractMeasure)).
"""
Statistics.std(μ::AbstractMeasure) = _sqrt_moment(var(μ))

"""
    cov(μ::AbstractMeasure)

The covariance matrix of the normalized law of `μ`, for array variates
(of their `vec`), see [`mean`](@ref Statistics.mean(::AbstractMeasure)).
Products and powers have known covariances only for scalar marginals.
"""
Statistics.cov(μ::AbstractMeasure) = _no_moment(μ, :covariance)

@noinline _no_moment(μ, what) =
    throw(ArgumentError("The $what of measures of type $(nameof(typeof(μ))) is not known analytically"))

@inline _sqrt_moment(v::Number) = sqrt(v)
@inline _sqrt_moment(v) = map(_sqrt_moment, v)

@inline _zero_moment(x::Number) = zero(float(x))
@inline _zero_moment(x) = map(_zero_moment, x)


Statistics.mean(::StdUniform) = 0.5
Statistics.var(::StdUniform) = 1 / 12
Statistics.mean(::StdNormal) = 0.0
Statistics.var(::StdNormal) = 1.0
Statistics.mean(::StdExponential) = 1.0
Statistics.var(::StdExponential) = 1.0
Statistics.mean(::StdLogistic) = 0.0
Statistics.var(::StdLogistic) = π^2 / 3

function Statistics.mean(t::TruncatedMeasure{<:StdNormal})
    first(_normal_truncated_moments(t.lo, t.hi))
end
function Statistics.var(t::TruncatedMeasure{<:StdNormal})
    last(_normal_truncated_moments(t.lo, t.hi))
end

# Moments of the conditional law, independent of the retained mass. Work
# at least in Float64: tail variances subtract terms of order one, even
# when the variance itself is small. Return the precision of the bounds.
function _normal_truncated_moments(lo, hi)
    a, b = promote(float(lo), float(hi))
    m, v = _normal_truncated_moments_impl(a * 1.0, b * 1.0)
    return oftype(a, m), oftype(a, v)
end

function _normal_truncated_moments_impl(a, b)
    a < b || return (oftype(a, NaN), oftype(a, NaN))
    if abs(a) > abs(b)
        m, v = _normal_truncated_moments_impl(-b, -a)
        return -m, v
    elseif isinf(a)
        return zero(a), one(a)
    end
    short = @_adignore (b - a) * max(abs((a + b) / 2), one(a)) <= one(a)
    if short
        return _normal_short_moments((a + b) / 2, (b - a) / 2)
    end
    mass = _std_normal_smf(a, b)
    p = exp(-a^2 / 2) * invsqrt2π / mass
    q = isinf(b) ? zero(p) : exp(-b^2 / 2) * invsqrt2π / mass
    m = p - q
    v = one(m) + (a - m) * p - (isinf(b) ? zero(q) : (b - m) * q)
    return m, v
end

# Integrate the density about the midpoint, on [-1, 1], without subtracting
# raw moments. Its coefficients obey n*c_n = -c*h*c_(n-1) - h^2*c_(n-2).
# With |c*h| <= 1/2 and h^2 <= 1/4, the omitted integrated absolute
# coefficients are bounded by Σ_{j+2k>20} (1/2)^j*(1/8)^k/(j!*k!*(j+2k+1))
# < 6e-17. The wider cutoff than for interval masses avoids cancellation
# of the centered second moment in the tail.
function _normal_short_moments(c, h)
    ch, h2 = c * h, h^2
    previous, coefficient = one(ch), -ch
    s0, s1, s2 = one(ch), coefficient / 3, one(ch) / 3
    for n in 2:20
        previous, coefficient = coefficient, (-ch * coefficient - h2 * previous) / n
        if iseven(n)
            s0 += coefficient / (n + 1)
            s2 += coefficient / (n + 3)
        else
            s1 += coefficient / (n + 2)
        end
    end
    m = s1 / s0
    return c + h * m, h2 * (s2 / s0 - m^2)
end

Statistics.mean(μ::Dirac) = μ.x
Statistics.var(μ::Dirac) = _zero_moment(μ.x)
function Statistics.cov(μ::Dirac{<:AbstractArray{<:Number}})
    n = length(μ.x)
    fill!(similar(μ.x, float(eltype(μ.x)), (n, n)), 0)
end

# Weights don't change the normalized law:
Statistics.mean(μ::AbstractWeightedMeasure) = mean(basemeasure(μ))
Statistics.var(μ::AbstractWeightedMeasure) = var(basemeasure(μ))
Statistics.cov(μ::AbstractWeightedMeasure) = cov(basemeasure(μ))

# The marginals of products are independent:
Statistics.mean(μ::AbstractProductMeasure) = map(mean, marginals(μ))
Statistics.var(μ::AbstractProductMeasure) = map(var, marginals(μ))
Statistics.cov(μ::AbstractProductMeasure) = _diagonal_cov(μ, var(μ))

@inline _diagonal_cov(::Any, v::AbstractArray{<:Number}) = Diagonal(vec(v))
@inline _diagonal_cov(μ, ::Any) = _no_moment(μ, :covariance)
