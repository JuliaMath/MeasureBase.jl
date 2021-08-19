export Affine, AffineTransform

struct AffineTransform{N,T}
    par::NamedTuple{N,T}
end

@inline Base.getproperty(d::AffineTransform, s::Symbol) = getfield(getfield(d, :par), s)

Base.propertynames(d::AffineTransform{N}) where {N} = N 

@inline Base.inv(f::AffineTransform{(:μ,:σ)}) = AffineTransform((μ = -(f.σ \ f.μ), ω = f.σ))
@inline Base.inv(f::AffineTransform{(:μ,:ω)}) = AffineTransform((μ = - f.ω * f.μ, σ = f.ω))

(f::AffineTransform{(:μ,:σ)})(x) = f.σ * x + f.μ

(f::AffineTransform{(:μ,:ω)})(x) = f.ω \ x + f.μ

###############################################################################

struct Affine{N,M,T} <: AbstractMeasure
    f::AffineTransform{N,T}
    parent::M
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

logdensity(d::Affine{(:μ,:σ)}, x) = logdensity(d.parent, d.σ \ (x - d.μ)) - log(d.σ)

logdensity(d::Affine{(:μ,:ω)}, x) = logdensity(d.parent, d.ω * (x - d.μ)) + log(d.ω)
