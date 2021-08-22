export Affine, AffineTransform

struct AffineTransform{N,T}
    par::NamedTuple{N,T}
end

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

###############################################################################

struct Affine{N,M,T} <: AbstractMeasure
    f::AffineTransform{N,T}
    parent::M
end

params(::Type{A}, nt::NamedTuple{C}) where {A<:Affine{N,M}} where {N,M,C} = tuple(setdiff(union(N, params(M)),C)...)


parent(d::Affine) = getfield(d, :parent)

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

basemeasure(d::Affine{N,L}) where {N, L<:Lebesgue} = d.parent

logdensity(d::Affine{N,L}, x) where {N,L<:Lebesgue} = logjac(getfield(d, :f))

# TODO: `log` doesn't work for the multivariate case, we need the log absolute determinant
logjac(::AffineTransform{(:μ,:σ)}) = -log(d.σ)
logjac(::AffineTransform{(:μ,:ω)}) = log(d.ω)
logjac(::AffineTransform{(:σ,)}) = -log(d.σ)
logjac(::AffineTransform{(:ω,)}) = log(d.ω)
logjac(::AffineTransform{(:μ,)}) = 0.0

function Base.rand(rng::Random.AbstractRNG, ::Type{T}, d::Affine) where {T}
    z = rand(rng, T, parent(d))
    f = getfield(d, :f)
    return f(z)
end
