export IntegerRange

abstract type AbstractDomain end

struct RealBounds{L,U} <: AbstractDomain
    lower :: L
    upper :: U
end

Base.in(x, b::RealBounds) = b.lower ≤ x ≤ b.upper


const ℝ = RealBounds(static(-Inf), static(Inf))

const ℝ₊ = RealBounds(static(0.0), static(Inf))

const 𝕀 = RealBounds(static(0.0), static(1.0))

Base.minimum(b::RealBounds) = b.lower
Base.maximum(b::RealBounds) = b.upper

Base.show(io::IO, ::typeof(ℝ)) = print(io, "ℝ")
Base.show(io::IO, ::typeof(ℝ₊)) = print(io, "ℝ₊")
Base.show(io::IO, ::typeof(𝕀)) = print(io, "𝕀")

testvalue(::typeof(ℝ)) = 0.0
testvalue(::typeof(ℝ₊)) = 1.0
testvalue(::typeof(𝕀)) = 0.5

struct IntegerBounds{L,U} <: AbstractDomain
    lower :: L
    upper :: U
end

Base.in(x, b::IntegerBounds) = b.lower ≤ x ≤ b.upper

const ℤ = IntegerBounds(static(-Inf), static(Inf))

Base.show(io::IO, ::typeof(ℤ)) = print(io, "ℤ")


###########################################################
# Integer ranges

struct IntegerRange{lo,hi} <: AbstractDomain end

Base.minimum(::IntegerRange{lo,hi}) where {lo,hi} = lo
Base.maximum(::IntegerRange{lo,hi}) where {lo,hi} = hi

Base.iterate(r::IntegerRange{lo,hi}) where {lo,hi} = iterate(lo:hi)

function Base.getindex(::Integers, r::AbstractUnitRange)
    IntegerRange{minimum(r),maximum(r)}()
end

function Base.show(io::IO, r::IntegerRange{lo,hi}) where {lo,hi}
    io = IOContext(io, :compact => true)
    print(io, "ℤ[", lo, ":", hi, "]")
end

testvalue(::IntegerRange{lo,hi}) where {lo,hi} = lo

###########################################################
# Real intervals

struct RealInterval{lo,hi} <: AbstractDomain end


struct Simplex{D} <: AbstractDomain
    dim::D # dimensionality as a manifold
end

projectto!(x, ::Simplex) = normalize!(x, 1)

struct Sphere{D} <: AbstractDomain
    dim::D # dimensionality as a manifold
end


projectto!(x, ::Sphere) = normalize!(x, 2)
struct ZeroSet{F, G} <: AbstractDomain
    f::F
    ∇f::G
end

function zeroset(::Simplex)
    f(x::AbstractVector{T}) where {T} = sum(x) - one(T)
    ∇f(x::AbstractVector{T}) where {T} = Fill(one(T), axes(x))
    ZeroSet(f, ∇f)
end

function zeroset(::Sphere)
    f(x::AbstractVector{T}) where {T} = (sum(xⱼ -> xⱼ^2, x) - one(T)) / 2
    ∇f(x::AbstractVector{T}) where {T} = x
    ZeroSet(f, ∇f)
end

struct LebesgueCodimOne{D,T,O} <: AbstractMeasure
    ndims ::D
    ortho ::O
end

function logdensityof(d::Density{L1, L2}, x) where {L1<:LebesgueCodimOne, L2<:LebesgueCodimOne}
    μ = d.μ
    ν = d.base
    μ.ndims == ν.ndims || return NaN
    rank([μ.ortho ν.ortho]) == 1 || return NaN
    return 0.0
end

struct LebesgueSimplex <: AbstractMeasure end

basemeasure(::LebesgueSimplex, x)