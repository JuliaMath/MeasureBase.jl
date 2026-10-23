# Moments of the normalized laws of measures, where they are known exactly.

"""
    mean(μ::AbstractMeasure)

The mean of the normalized law of `μ`.

Moments are known exactly for wrapped distributions, standard measures,
point masses, products and powers (in the layout of their variates),
weighted and normalized measures of measures with known moments and, with
AffineMaps loaded, affine pushforwards of those. Other measures throw an
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
