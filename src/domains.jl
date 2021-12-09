abstract type AbstractDomain end

abstract type RealDomain <: AbstractDomain end

# TODO: Use IntervalSets
struct RealNumbers <: RealDomain end

const ℝ = RealNumbers()

Base.minimum(::RealNumbers) = static(-Inf)
Base.maximum(::RealNumbers) = static(Inf)

Base.in(x, ::RealNumbers) = isreal(x)

Base.show(io::IO, ::typeof(ℝ)) = print(io, "ℝ")

struct BoundedReals{L,U} <: RealDomain
    lower::L
    upper::U
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
    lower::L
    upper::U
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
# ZeroSet

struct ZeroSet{F, G} <: AbstractDomain
    f::F
    ∇f::G
end

# Based on some quick tests, but may need some adjustment
Base.in(x::AbstractArray{T}, z::ZeroSet) where {T} = abs(z.f(x)) < ldexp(eps(float(T)), 6)


###########################################################
# CodimOne

abstract type CodimOne <: AbstractDomain end

function tangentat(a::CodimOne, b::CodimOne, x::AbstractArray{T}; tol=ldexp(eps(float(T)), 6)) where {T}
    # Sometimes you get lucky
    a == b && return true

    # Get the normal vectors
    g1 = a.∇f(x)
    g2 = b.∇f(x)
    
    # See if one is a multiple of the other
    one(T) - Statistics.corm(g1, zero(T), g2, zero(T)) < tol
end

function zeroset(::CodimOne)::ZeroSet end

###########################################################
# Simplex

struct Simplex <: CodimOne end

function zeroset(::Simplex) 
    f(x::AbstractArray{T}) where {T} = sum(x) - one(T)
    ∇f(x::AbstractArray{T}) where {T} = Fill(one(T), size(x))
    ZeroSet(f, ∇f)
end

function Base.in(x::AbstractArray{T}, ::Simplex) where {T} 
    x .≥ zero(eltype(x)) || return false
    return x ∈ zeroset(Simplex())
end

projectto!(x, ::Simplex) = normalize!(x, 1)

###########################################################
# Sphere

struct Sphere <: CodimOne end

function zeroset(::Sphere) 
    f(x::AbstractArray{T}) where {T} = sum(xⱼ -> xⱼ^2, x) - one(T)
    ∇f(x::AbstractArray{T}) where {T} = x
    ZeroSet(f, ∇f)
end

function Base.in(x::AbstractArray{T}, ::Sphere) where {T} 
    return x ∈ zeroset(Sphere())
end

projectto!(x, ::Sphere) = normalize!(x, 2)

