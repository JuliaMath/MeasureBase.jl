# This file is a part of DistributionMeasures.jl, licensed under the MIT License (MIT).

"""
    const StandardNormal{N} = StandardDist{Normal,N}

The univariate standard normal distribution.
"""
const StandardNormal{N} = StandardDist{Normal,N}
export StandardNormal

Distributions.Normal(d::StandardDist{Normal,0}) = Distributions.Normal()
Base.convert(::Type{Distributions.Normal}, d::StandardDist{Normal,1}) = Distributions.Normal(d)

Distributions.MvNormal(d::StandardDist{Normal,1}) = MvNormal(PDMats.ScalMat(length(d), 1))
Base.convert(::Type{Distributions.MvNormal}, d::StandardDist{Normal,1}) = Distributions.MvNormal(d)

Base.minimum(d::StandardDist{Normal,0}) = -Inf
Base.maximum(d::StandardDist{Normal,0}) = +Inf

Distributions.insupport(d::StandardDist{Normal,0}, x::Real) = !isnan(x)

Distributions.location(d::StandardDist{Normal,0}) = Statistics.mean(d)
Distributions.scale(d::StandardDist{Normal,0}) = Statistics.var(d)

Statistics.mean(d::StandardDist{Normal,0}) = 0
Statistics.mean(d::StandardDist{Normal,N}) where N = Zeros{Int}(size(d)...)

StatsBase.median(d::StandardDist{Normal}) = Statistics.mean(d)
StatsBase.mode(d::StandardDist{Normal}) = Statistics.mean(d)

StatsBase.modes(d::StandardDist{Normal,0}) = Zeros{Int}(1)

Statistics.var(d::StandardDist{Normal,0}) = 1
Statistics.var(d::StandardDist{Normal,N}) where N = Ones{Int}(size(d)...)

StatsBase.std(d::StandardDist{Normal,0}) = 1
StatsBase.std(d::StandardDist{Normal,N}) where N = Ones{Int}(size(d)...)

StatsBase.skewness(d::StandardDist{Normal,0}) = 0
StatsBase.kurtosis(d::StandardDist{Normal,0}) = 0

StatsBase.entropy(d::StandardDist{Normal,0}) = muladd(log2π, 1/2, 1/2)

Distributions.logpdf(d::StandardDist{Normal,0}, x::U) where {U<:Real} = muladd(abs2(x), -U(1)/U(2), -log2π/U(2))
Distributions.pdf(d::StandardDist{Normal,0}, x::U) where {U<:Real} = invsqrt2π * exp(-abs2(x)/U(2))

@inline Distributions.gradlogpdf(d::StandardDist{Normal,0}, x::Real) = -x

@inline Distributions.logcdf(d::StandardDist{Normal,0}, x::Real) = StatsFuns.normlogcdf(x)
@inline Distributions.cdf(d::StandardDist{Normal,0}, x::Real) = StatsFuns.normcdf(x)
@inline Distributions.logccdf(d::StandardDist{Normal,0}, x::Real) = StatsFuns.normlogccdf(x)
@inline Distributions.ccdf(d::StandardDist{Normal,0}, x::Real) = StatsFuns.normccdf(x)
@inline Distributions.quantile(d::StandardDist{Normal,0}, p::Real) = StatsFuns.norminvcdf(p)
@inline Distributions.cquantile(d::StandardDist{Normal,0}, p::Real) = StatsFuns.norminvccdf(p)
@inline Distributions.invlogcdf(d::StandardDist{Normal,0}, p::Real) = StatsFuns.norminvlogcdf(p)
@inline Distributions.invlogccdf(d::StandardDist{Normal,0}, p::Real) = StatsFuns.norminvlogccdf(p)

Base.rand(rng::AbstractRNG, d::StandardDist{Normal,0}) = randn(rng)
Base.rand(rng::AbstractRNG, d::StandardDist{Normal,N}) where N = randn(rng, size(d)...)
Random.rand!(rng::AbstractRNG, d::StandardDist{Normal,N}, x::AbstractArray{<:Real,N}) where N = Random.randn!(rng, x)

Distributions.invcov(d::StandardDist{Normal,1}) = Distributions.cov(d)
Distributions.logdetcov(d::StandardDist{Normal,1}) = 0


function Distributions.sqmahal(d::StandardDist{Normal,N}, x::AbstractArray{<:Real,N}) where N
    dot(x, checked_arg(d, x))
end

function Distributions. sqmahal!(r::AbstractVector, d::StandardDist{Normal,N}, x::AbstractMatrix) where N
    x_cols = eachcol(checked_arg(d, first(eachcol(x))))
    r .= dot.(x_cols, x_cols)
end
