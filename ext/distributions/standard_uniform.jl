# This file is a part of DistributionMeasures.jl, licensed under the MIT License (MIT).

"""
    const StandardUniform{N} = StandardDist{Uniform,N}

The univariate standard uniform distribution.
"""
const StandardUniform{N} = StandardDist{Uniform,N}
export StandardUniform

Distributions.Uniform(d::StandardDist{Uniform,0}) = Distributions.Uniform()
Base.convert(::Type{Distributions.Uniform}, d::StandardDist{Uniform,1}) = Distributions.Uniform(d)

Base.minimum(::StandardDist{Uniform,0}) = 0
Base.maximum(::StandardDist{Uniform,0}) = 1

Distributions.location(::StandardDist{Uniform,0}) = 0
Distributions.scale(::StandardDist{Uniform,0}) = 1

Statistics.mean(d::StandardDist{Uniform,0}) = 1//2
StatsBase.median(d::StandardDist{Uniform,0}) = Statistics.mean(d)
StatsBase.mode(d::StandardDist{Uniform,0}) = Statistics.mean(d)
StatsBase.modes(d::StandardDist{Uniform,0}) = Zeros{Int}(0)
StatsBase.modes(d::StandardDist{Uniform,N}) where N = Fill(Zeros{Int}(size(d)))

Statistics.var(d::StandardDist{Uniform,0}) = 1//12
StatsBase.std(d::StandardDist{Uniform,0}) = sqrt(Statistics.var(d))
StatsBase.skewness(d::StandardDist{Uniform,0}) = 0
StatsBase.kurtosis(d::StandardDist{Uniform,0}) = -6//5

StatsBase.entropy(d::StandardDist{Uniform,0}) = 0


function Distributions.logpdf(d::StandardDist{Uniform,0}, x::U) where {U<:Real}
    ifelse(Distributions.insupport(d, x), U(0), U(-Inf))
end

function Distributions.pdf(d::StandardDist{Uniform,0}, x::U) where {U<:Real}
    ifelse(Distributions.insupport(d, x), one(U), zero(U))
end


Distributions.logcdf(d::StandardDist{Uniform,0}, x::U) where {U<:Real} = log(Distributions.cdf(d, x))

function Distributions.cdf(d::StandardDist{Uniform,0}, x::U) where {U<:Real}
    ifelse(x < zero(U), zero(U), ifelse(x < one(U), x, one(U)))
end

Distributions.logccdf(d::StandardDist{Uniform,0}, x::U) where {U<:Real} = log(Distributions.ccdf(d, x))

Distributions.ccdf(d::StandardDist{Uniform,0}, x::U) where {U<:Real} = one(x) - Distributions.cdf(d, x)


function Distributions.quantile(d::StandardDist{Uniform,0}, p::U) where {U<:Real}
   convert(float(U), p)
end

function Distributions.cquantile(d::StandardDist{Uniform,0}, p::U) where {U<:Real}
    y = Distributions.quantile(d, p)
    one(y) - y
end


Distributions.mgf(d::StandardDist{Uniform,0}, t::Real) = Distributions.mgf(nonstddist(d), t)
Distributions.cf(d::StandardDist{Uniform,0}, t::Real) = Distributions.cf(nonstddist(d), t)

Distributions.gradlogpdf(d::StandardDist{Uniform,0}, x::Real) = zero(x)

function Distributions.gradlogpdf(d::StandardDist{Uniform,N}, x::AbstractArray{<:Real,N}) where N
    zero(checked_arg(d, x))
end

Base.rand(rng::AbstractRNG, d::StandardDist{Uniform,0}) = rand(rng)
Base.rand(rng::AbstractRNG, d::StandardDist{Uniform,N}) where N = rand(rng, size(d)...)
Random.rand!(rng::AbstractRNG, d::StandardDist{Uniform,N}, x::AbstractArray{<:Real,N}) where N = rand!(rng, x)
