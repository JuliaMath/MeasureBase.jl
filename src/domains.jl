"""Abstract supertype for all measure domains."""
abstract type AbstractDomain end

"""
    @domain(name, T)

Defines a new singleton struct `T`, and a value `name` for building values of
that type.

# Examples

For instance, `@domain ℝ RealNumbers` is equivalent to
```julia
struct RealNumbers <: AbstractDomain end
ℝ = RealNumbers()
Base.show(io::IO, ::RealNumbers) = print(io, "ℝ")
```
"""
macro domain(name, T)
    sname = String(name)

    name = esc(name)
    quote
        struct $T <: AbstractDomain end
        const $name = $T()
        Pretty.tile(::$T) = Pretty.literal($sname)
    end
end

@domain ℝ RealNumbers

@domain ℝ₊ PositiveReals

@domain 𝕀 UnitInterval

@domain ℤ Integers

###########################################################
# Integer ranges

"""
    IntegerRange{lo,hi}

Domain containing all the integers between `lo` and `hi` (inclusive).
"""
struct IntegerRange{lo,hi} <: AbstractDomain end

Base.minimum(::IntegerRange{lo,hi}) where {lo,hi} = lo
Base.maximum(::IntegerRange{lo,hi}) where {lo,hi} = hi

Base.iterate(::IntegerRange{lo,hi}) where {lo,hi} = iterate(lo:hi)

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

"""
    RealInterval{lo,hi}

Domain containing all the real numbers between `lo` and `hi` (inclusive).

!!! warning "Missing methods"
    Methods for this domain are not yet implemented. See [this issue](https://github.com/cscherrer/MeasureBase.jl/issues/16) for the discussion.
"""
struct RealInterval{lo,hi} <: AbstractDomain end
