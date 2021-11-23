export IntegerRange

abstract type AbstractDomain end

"""
    @domain(name, T)

Defines a new singleton struct `T`, and a value `name` for building values of
that type.

For example, `@domain ℝ RealNumbers` is equivalent to

    struct RealNumbers <: AbstractDomain end

    export ℝ

    ℝ = RealNumbers()

    Base.show(io::IO, ::RealNumbers) = print(io, "ℝ")
"""
macro domain(name, T)
    sname = String(name)

    name = esc(name)
    quote
        struct $T <: AbstractDomain end
        export $name
        const $name = $T()
        Base.show(io::IO, ::$T) = Pretty.literal($sname)
    end
end

@domain ℝ RealNumbers

@domain ℝ₊ PositiveReals

@domain 𝕀 UnitInterval

@domain ℤ Integers

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


