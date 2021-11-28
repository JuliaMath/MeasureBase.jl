abstract type AbstractDomain end

abstract type RealDomain <: AbstractDomain end

struct RealNumbers <: RealDomain end

const ℝ = RealNumbers()

Base.minimum(::RealNumbers) = static(-Inf)
Base.maximum(::RealNumbers) = static(Inf)

Base.in(x, ::RealNumbers) = isreal(x)

Base.show(io::IO, ::typeof(ℝ)) = print(io, "ℝ")

struct BoundedReals{L,U} <: RealDomain
    lower :: L
    upper :: U
end

Base.in(x, b::BoundedReals) = b.lower ≤ x ≤ b.upper


export ℝ, ℝ₊, 𝕀, ℤ


const ℝ₊ = BoundedReals(static(0.0), static(Inf))
const 𝕀 = BoundedReals(static(0.0), static(1.0))

Base.minimum(b::BoundedReals) = b.lower
Base.maximum(b::BoundedReals) = b.upper


Base.show(io::IO, ::typeof(ℝ₊)) = print(io, "ℝ₊")
Base.show(io::IO, ::typeof(𝕀)) = print(io, "𝕀")

testvalue(::typeof(ℝ)) = 0.0
testvalue(::typeof(ℝ₊)) = 1.0
testvalue(::typeof(𝕀)) = 0.5

abstract type IntegerDomain <: AbstractDomain end

struct IntegerNumbers <: IntegerDomain end

Base.in(x, ::IntegerNumbers) = isinteger(x)

const ℤ = IntegerNumbers()

Base.show(io::IO, ::typeof(ℤ)) = print(io, "ℤ")

Base.minimum(::IntegerNumbers) = static(-Inf)
Base.maximum(::IntegerNumbers) = static(Inf)
struct BoundedInts{L,U} <: IntegerDomain
    lower :: L
    upper :: U
end

Base.in(x, b::BoundedInts) = x ∈ ℤ && b.lower ≤ x ≤ b.upper

Base.minimum(b::BoundedInts) = b.lower
Base.maximum(b::BoundedInts) = b.upper

function Base.show(io::IO, b::BoundedInts)
    io = IOContext(io, :compact => true)
    print(io, "ℤ[", b.lower, ":", b.upper, "]")
end

testvalue(b::BoundedInts) = min(b.lower, 0)

function Base.getindex(::typeof(ℤ), r::AbstractUnitRange)
    BoundedInts(extrema(r)...)
end





###########################################################
# Simplex

# struct Simplex{D} <: AbstractDomain
#     dim::D # dimensionality as a manifold
# end

# projectto!(x, ::Simplex) = normalize!(x, 1)

# struct Sphere{D} <: AbstractDomain
#     dim::D # dimensionality as a manifold
# end


# projectto!(x, ::Sphere) = normalize!(x, 2)
# struct ZeroSet{F, G} <: AbstractDomain
#     f::F
#     ∇f::G
# end

# function zeroset(::Simplex)
#     f(x::AbstractVector{T}) where {T} = sum(x) - one(T)
#     ∇f(x::AbstractVector{T}) where {T} = Fill(one(T), axes(x))
#     ZeroSet(f, ∇f)
# end

# function zeroset(::Sphere)
#     f(x::AbstractVector{T}) where {T} = (sum(xⱼ -> xⱼ^2, x) - one(T)) / 2
#     ∇f(x::AbstractVector{T}) where {T} = x
#     ZeroSet(f, ∇f)
# end

# struct LebesgueCodimOne{D,T,O} <: AbstractMeasure
#     ndims ::D
#     ortho ::O
# end

# function logdensityof(d::Density{L1, L2}, x) where {L1<:LebesgueCodimOne, L2<:LebesgueCodimOne}
#     μ = d.μ
#     ν = d.base
#     μ.ndims == ν.ndims || return NaN
#     rank([μ.ortho ν.ortho]) == 1 || return NaN
#     return 0.0
# end

# struct LebesgueSimplex <: AbstractMeasure end

# basemeasure(::LebesgueSimplex, x)