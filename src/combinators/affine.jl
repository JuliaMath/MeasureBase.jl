export Affine, AffineTransform

struct AffineTransform{N,T}
    par::NamedTuple{N,T}
end

params(f::AffineTransform) = getfield(f, :par)

@inline Base.getproperty(d::AffineTransform, s::Symbol) = getfield(getfield(d, :par), s)

Base.propertynames(d::AffineTransform{N}) where {N} = N 

@inline Base.inv(f::AffineTransform{(:μ,:σ)}) = AffineTransform((μ = -(f.σ \ f.μ), ω = f.σ))
@inline Base.inv(f::AffineTransform{(:μ,:ω)}) = AffineTransform((μ = - f.ω * f.μ, σ = f.ω))
@inline Base.inv(f::AffineTransform{(:σ,)}) = AffineTransform((ω = f.σ,))
@inline Base.inv(f::AffineTransform{(:ω,)}) = AffineTransform((σ = f.ω,))
@inline Base.inv(f::AffineTransform{(:μ,)}) = AffineTransform((μ = -f.μ,))

(f::AffineTransform{(:μ,)})(x) = x + f.μ
(f::AffineTransform{(:σ,)})(x) = f.σ * x
(f::AffineTransform{(:ω,)})(x) = f.ω \ x
(f::AffineTransform{(:μ,:σ)})(x) = f.σ * x + f.μ
(f::AffineTransform{(:μ,:ω)})(x) = f.ω \ x + f.μ

# TODO: `log` doesn't work for the multivariate case, we need the log absolute determinant
logjac(f::AffineTransform{(:μ,:σ)}) = log(f.σ)
logjac(f::AffineTransform{(:μ,:ω)}) = -log(f.ω)
logjac(f::AffineTransform{(:σ,)}) = log(f.σ)
logjac(f::AffineTransform{(:ω,)}) = -log(f.ω)
logjac(f::AffineTransform{(:μ,)}) = 0.0

###############################################################################

struct Affine{N,M,T} <: AbstractMeasure
    f::AffineTransform{N,T}
    parent::M

    function Affine(f::AffineTransform, parent::WeightedMeasure)
        WeightedMeasure(parent.logweight, Affine(f, parent.base))
    end

    Affine(f::AffineTransform{N,T}, parent::M) where {N,M,T} = new{N,M,T}(f, parent)
end

parent(d::Affine) = getfield(d, :parent)

function params(μ::Affine)
    nt1 = getfield(getfield(μ, :f), :par)
    nt2 = params(parent(μ))
    return merge(nt1, nt2)
end

function paramnames(::Type{A}) where {N,M, A<:Affine{N,M}}
    tuple(union(N, paramnames(M))...)
end

Affine(nt::NamedTuple, μ::AbstractMeasure) = Affine(AffineTransform(nt), μ)

Affine(nt::NamedTuple) = μ -> Affine(nt, μ)

Base.propertynames(d::Affine{N}) where {N} = N ∪ (:parent,)

@inline function Base.getproperty(d::Affine, s::Symbol) 
    if s === :parent
        return getfield(d, :parent)
    else
        return getproperty(getfield(d, :f), s)
    end
end

# Note: We could also write
# logdensity(d::Affine, x) = logdensity(inv(getfield(d, :f)), x)

logdensity(d::Affine{(:μ,:σ)}, x) = logdensity(d.parent, d.σ \ (x - d.μ))
logdensity(d::Affine{(:μ,:ω)}, x) = logdensity(d.parent, d.ω * (x - d.μ))
logdensity(d::Affine{(:σ,)}, x) = logdensity(d.parent, d.σ \ x)
logdensity(d::Affine{(:ω,)}, x) = logdensity(d.parent, d.ω * x)
logdensity(d::Affine{(:μ,)}, x) = logdensity(d.parent, x - d.μ) 

basemeasure(d::Affine) = Affine(getfield(d, :f), basemeasure(d.parent))

# We can't do this until we know we're working with Lebesgue measure, since for
# example it wouldn't make sense to apply a log-Jacobian to a point measure
basemeasure(d::Affine{N,L}) where {N, L<:Lebesgue} = WeightedMeasure(-logjac(d), d.parent)

logjac(d::Affine) = logjac(getfield(d, :f))


function Base.rand(rng::Random.AbstractRNG, ::Type{T}, d::Affine) where {T}
    z = rand(rng, T, parent(d))
    f = getfield(d, :f)
    return f(z)
end
