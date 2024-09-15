# This file is a part of DistributionMeasures.jl, licensed under the MIT License (MIT).

"""
    struct StandardDist{D<:Distribution{Univariate,Continuous},N} <: Distributions.Distribution{ArrayLikeVariate{N},Continuous}

Represents `D()` or a product distribution of `D()` in a dispatchable fashion.

Constructor:
```
    StandardDist{Uniform}(size...)
    StandardDist{Normal}(size...)
```
"""
struct StandardDist{D<:Distribution{Univariate,Continuous},N,U<:Integer} <: Distributions.Distribution{ArrayLikeVariate{N},Continuous}
    _size::NTuple{N,U}
end
export StandardDist

StandardDist{D}() where {D<:Distribution{Univariate,Continuous}} = StandardDist{D,0,StaticInt{1}}(())
StandardDist{D}(dims::Vararg{U,N}) where {D<:Distribution{Univariate,Continuous},N,U<:Integer} = StandardDist{D,N,U}((dims...,))


const StandardUnivariateDist{D<:Distribution{Univariate,Continuous},U<:Integer} = StandardDist{D,0,U}
const StandardMultivariteDist{D<:Distribution{Univariate,Continuous},U<:Integer} = StandardDist{D,1,U}


function Base.show(io::IO, d::StandardDist{D}) where D
    print(io, nameof(typeof(d)), "{", D, "}")
    show(io, d._size)
end


@inline MeasureBase.transport_def(::MU, μ::MU, x) where {MU<:StandardDist{<:Any,0}} = x

for (A, B) in [
    (Uniform, StdUniform),
    (Exponential, StdExponential),
    (Logistic, StdLogistic),
    (Normal, StdNormal)
]
    @eval begin
        @inline MeasureBase.transport_origin(d::StandardDist{$A,0}) = $B()
        @inline MeasureBase.transport_origin(d::StandardDist{$A,N}) where N = $B()^size(d)
    end
end

@inline MeasureBase.to_origin(ν::StandardDist, y) = y
@inline MeasureBase.from_origin(ν::StandardDist, x) = x


@inline nonstddist(::StandardDist{D,0}) where D = D(Distributions.params(D())...)
@inline function nonstddist(d::StandardDist{D,N}) where {D,N}
    nonstd0 = nonstddist(StandardDist{D}())
    reshape(Distributions.product_distribution(fill(nonstd0, length(d))), size(d))
end


(::Type{D}, d::StandardDist{D,0}) where {D<:Distribution{Univariate,Continuous}} = nonstddist(d)

# TODO: Replace `fill` by `FillArrays.Fill` once Distributions fully supports this:
(::Type{Distributions.Product})(d::StandardDist{D,1}) where D = Distributions.Product(fill(StandardDist{D}(), length(d)))

Base.convert(::Type{D}, d::StandardDist{D,0}) where {D<:Distribution{Univariate,Continuous}} = D(d)
Base.convert(::Type{Distributions.Product}, d::StandardDist{D,1}) where D = Distributions.Product(d)



@inline Base.size(d::StandardDist) = d._size
@inline Base.length(d::StandardDist) = prod(size(d))

Base.eltype(::Type{StandardDist{D,N}}) where {D,N} = Float64

@inline Distributions.partype(d::StandardDist{D}) where D = Float64

@inline StatsBase.params(d::StandardDist) = ()

for f in (
    :(Base.minimum),
    :(Base.maximum),
    :(Statistics.mean),
    :(Statistics.median),
    :(StatsBase.mode),
    :(Statistics.var),
    :(Statistics.std),
    :(StatsBase.skewness),
    :(StatsBase.kurtosis),
    :(Distributions.location),    
    :(Distributions.scale),    
)
    @eval begin
        ($f)(d::StandardDist{D,0}) where D = ($f)(nonstddist(d))
        ($f)(d::StandardDist{D,N}) where {D,N} = Fill(($f)(StandardDist{D}()), size(d)...)
    end
end

StatsBase.modes(d::StandardDist) = [StatsBase.mode(d)]

# ToDo: Define cov for N!=1?
Statistics.cov(d::StandardDist{D,1}) where D = Diagonal(Statistics.var(d))
Distributions.invcov(d::StandardDist{D,1}) where D = Diagonal(Fill(inv(Statistics.var(StandardDist{D}())), length(d)))
Distributions.logdetcov(d::StandardDist{D,1}) where D = length(d) * log(Statistics.var(StandardDist{D}()))

StatsBase.entropy(d::StandardDist{D,0}) where D = StatsBase.entropy(nonstddist(d))
StatsBase.entropy(d::StandardDist{D,N}) where {D,N} = length(d) * StatsBase.entropy(StandardDist{D}())


Distributions.insupport(d::StandardDist{D,0}, x::Real) where D = Distributions.insupport(nonstddist(d), x)

function Distributions.insupport(d::StandardDist{D,N}, x::AbstractArray{<:Real,N}) where {D,N}
    all(Base.Fix1(Distributions.insupport, StandardDist{D}()), checked_arg(d, x))
end


@inline Distributions.logpdf(d::StandardDist{D,0}, x::U) where {D,U} = Distributions.logpdf(nonstddist(d), x)

function Distributions.logpdf(d::StandardDist{D,N}, x::AbstractArray{<:Real,N}) where {D,N}
    Distributions._logpdf(d, checked_arg(d, x))
end

function Distributions._logpdf(::StandardDist{D,1}, x::AbstractArray{<:Real,1}) where D
    sum(Base.Fix1(Distributions.logpdf, StandardDist{D}()), x)
end

function Distributions._logpdf(::StandardDist{D,2}, x::AbstractArray{<:Real,2}) where D
    sum(Base.Fix1(Distributions.logpdf, StandardDist{D}()), x)
end

function Distributions._logpdf(::StandardDist{D,N}, x::AbstractArray{<:Real,N}) where {D,N}
    sum(Base.Fix1(Distributions.logpdf, StandardDist{D}()), x)
end



Distributions.gradlogpdf(d::StandardDist{D,0}, x::Real) where D = Distributions.gradlogpdf(nonstddist(d), x)

function Distributions.gradlogpdf(d::StandardDist{D,N}, x::AbstractArray{<:Real,N}) where {D,N}
    Distributions.gradlogpdf.(StandardDist{D}(), checked_arg(d, x))
end


#@inline Distributions.pdf(d::StandardDist{D,0}, x::U) where {D,U} = pdf(nonstddist(d), x)

function Distributions.pdf(d::StandardDist{D,1}, x::AbstractVector{U}) where {D,U<:Real}
    Distributions._pdf(d, checked_arg(d, x))
end

function Distributions._pdf(d::StandardDist{D,1}, x::AbstractVector{U}) where {D,U<:Real}
    exp(Distributions._logpdf(d, x))
end

function Distributions.pdf(d::StandardDist{D,2}, x::AbstractMatrix{U}) where {D,U<:Real}
    Distributions._pdf(d, checked_arg(d, x))
end

function Distributions._pdf(d::StandardDist{D,2}, x::AbstractMatrix{U}) where {D,U<:Real}
    exp(Distributions._logpdf(d, x))
end

function Distributions.pdf(d::StandardDist{D,N}, x::AbstractArray{U,N}) where {D,N,U<:Real}
    Distributions._pdf(d, checked_arg(d, x))
end

function Distributions._pdf(d::StandardDist{D,N}, x::AbstractArray{U,N}) where {D,N,U<:Real}
    exp(Distributions._logpdf(d, x))
end


for f in (
    :(Distributions.logcdf),
    :(Distributions.cdf),
    :(Distributions.logccdf),
    :(Distributions.ccdf),
    :(Distributions.quantile),
    :(Distributions.cquantile),
    :(Distributions.invlogcdf),
    :(Distributions.invlogccdf),
    :(Distributions.mgf),
    :(Distributions.cf),
)
    @eval begin
        @inline ($f)(d::StandardDist, x::Real) = ($f)(nonstddist(d), x)
    end
end


Base.rand(rng::AbstractRNG, d::StandardDist{D,0}) where D = rand(rng, nonstddist(d))
Random.rand!(rng::AbstractRNG, d::StandardDist{D,0}, x::AbstractArray{<:Real,0}) where D = (x[] = rand(rng, d); return x)
Random.rand!(rng::AbstractRNG, d::StandardDist{D,N}, x::AbstractArray{<:Real,N}) where {D,N} = rand!(rng, StandardDist{D}(), x)


Distributions.truncated(d::StandardDist{D,0}, l::Real, u::Real) where D = Distributions.truncated(nonstddist(d), l, u)

Distributions.product_distribution(dists::AbstractVector{<:StandardDist{D,0}}) where D = StandardDist{D}(size(dists)...)
Distributions.product_distribution(dists::AbstractArray{<:StandardDist{D,0}}) where D = StandardDist{D}(size(dists)...)
