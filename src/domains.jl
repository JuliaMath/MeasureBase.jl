# Custom abstract set type. Design reserve to be able to switch to
#`Base.AbstractSet` or another set type hierarchy in the future:
"""
    MeasureBase.MAbstractSet

Abstract type for some measurable sets.

Not every measurable set needs to be a a subtype of
`MeasureBase.MAbstractSet`.

See also [`MeasureBase.SetLike`](@ref).
"""
abstract type MAbstractSet end

"""
    const MeasureBase.SetLike = Union{MeaureBase.MAbstractSet, Base.AbstractSet, IntervalSets.Domain}

Any kind of (measurable) set.

There needs to be an implicit sigma-algebra for subtypes of
`MeasureBase.SetLike` to make them useable for measures. This can't easily be
imposed via type constraints, though, to is is by-contract.
"""
const SetLike = Union{MeasureBase.MAbstractSet,Base.AbstractSet,IntervalSets.Domain}

"""
    mdomain(m)::MeasureBase.SetLikeIn

Return the domain, i.e. the measurable set, of the measure `m`.

The measure must allow of evaluating densities and the like over the whole
domain, even if the support of the measure is only a subset of the domain.

May return [`MeasureBase.ImplicitDomain`](@ref) if the domain cannot be
computed (efficiently).
"""
function mdomain end

@inline mdomain(m) = ImplicitDomain(m)

"""
    maybe_in(x, s)

Test if `x` may be a member of `s`.

Defaults to `in(x, s)`, but may be specialized for certain types of `s`,
e.g. for `s::MeasureBase.ImplicitDomain`.
"""
function maybe_in end

maybe_in(x, s) = in(x, s)

"""
    struct MeasureBase.ImplicitDomain{M} <: MeasureBase.MAbstractSet

Represents the domain (i.e. the measurable set) of a measure `m::M`.

Constructors:

```
MeasureBase.ImplicitDomain(m)
```

Fields:

* `m::M`: The measure.

For many pushforward measures and similar, the measureable space can not be
computed efficiently or at all. In such cases, [`mdomain(m)`](@ref) should
return `ImplicitDomain(m)`.

Does not support `Base.in(x, s::MeasureBase.ImplicitDomain)`, and
`MeasureBase.maybe_in(x, s::MeasureBase.ImplicitDomain)` always return `true`
(unless specialized for the measure type).
"""
struct ImplicitDomain{M} <: MAbstractSet
    m::M
end

@inline Base.union(s::ImplicitDomain, others::ImplicitDomain...) =
    ImplicitDomain(+(s.m, map(x -> x.m, others)...))

function Base.in(@nospecialize(x), ::ImplicitDomain)
    throw(ArgumentError("Cannot test if a value lies withing an implicit domain."))
end

maybe_in(@nospecialize(x), ::ImplicitDomain) = true

function Base.isempty(::ImplicitDomain)
    throw(ArgumentError("Can't test if an ImplicitDomain is empty"))
end

"""
    RealInterval() isa MeasureBase.MAbstractSet

The real numbers.
"""
struct RealNumbers <: MAbstractSet end

@inline Base.in(x::Real, ::RealNumbers) = true
@inline Base.in(x, ::RealNumbers) = isreal(x)

@inline Base.isempty(::RealNumbers) = false

@inline Base.union(s::RealNumbers, ::RealNumbers...) = s

@inline Base.minimum(::RealNumbers) = static(-Inf)
@inline Base.maximum(::RealNumbers) = static(Inf)

"""
    const MeasureBase.ℝ = RealNumbers()

The real numbers, see [`MeasureBase.RealNumbers`](@ref).
"""
const ℝ = RealNumbers()

Base.show(io::IO, ::MIME"text/plain", ::RealNumbers) = print(io, "MeasureBase.ℝ")

"""
    MeasureBase.Integers() isa MeasureBase.MAbstractSet
"""
struct Integers <: MAbstractSet end

@inline Base.in(x::Integer, ::Integers) = true
@inline Base.in(x, ::Integers) = isinteger(x)

@inline Base.isempty(::Integers) = false

@inline Base.union(s::Integers, ::Integers...) = s

# # This could get tricky with mixed-precision code. Probably needs some
# # special AbstractInteger infinity type (but custom AbstractInteger types
# # may cause a lot of method invalidations, which is why Static.StaticInteger
# # is not an AbstractInteger).
# @inline Base.minimum(::RealNumbers) = static(typemax(Int64))
# @inline Base.maximum(::RealNumbers) = static(typemin(Int64))

"""
    const ℤ = Integers()

The integers, see [`MeasureBase.Integers`](@ref).
"""
const ℤ = Integers()

Base.show(io::IO, ::MIME"text/plain", ::Integers) = print(io, "MeasureBase.ℤ")

"""
    struct MeasureBase.AbstractCartSetProd <: MAbstractSet

Supertype for cartesian products of sets.
"""
abstract type AbstractCartSetProd <: MAbstractSet end

"""
    struct SetCardProd <: AbstractCartSetProd

A cartesian product over a collection of sets.

Constructor:

```julia
prodset = SetCardProd(sets)
```

`sets` may be a `Tuple`, `NamedTuple` or `AbstractArray` of sets/domains.
"""
struct SetCardProd{S<:Union{Tuple,NamedTuple,AbstractArray}} <: AbstractCartSetProd
    sets::S
end

@inline Base.in(x::Tuple{Vararg{Any,N}}, s::SetCardProd{<:Tuple{Vararg{Any,N}}}) where {N} =
    prod(map(in, x, s.sets))::Bool
@inline Base.in(x::NamedTuple{names}, s::SetCardProd{<:NamedTuple{names}}) where {names} =
    prod(map(in, values(x), values(s.sets)))::Bool
# ToDo: Allow this?
# Base.in(x::AbstractVector, s::SetCardProd{<:Tuple}) = all(in.(x,s.sets))::Bool
function Base.in(
    x::AbstractArray{<:Any,N},
    s::SetCardProd{<:AbstractArray{<:Any,N}},
) where {N}
    all(in.(x, s.sets))::Bool
end

@inline Base.isempty(s::SetCardProd) = all(!isempty, s.sets)

@inline function Base.union(
    s::SetCardProd{<:Tuple{Vararg{Any,N}}},
    others::SetCardProd{<:Tuple{Vararg{Any,N}}}...,
) where {N}
    SetCardProd(map(union, s.sets, map(componentsets, others)...))
end

@inline function Base.union(
    s::SetCardProd{<:NamedTuple{names}},
    others::SetCardProd{<:NamedTuple{names}}...,
) where {names}
    SetCardProd(map(union, s.sets, map(componentsets, others)...))
end

function Base.union(
    s::SetCardProd{<:AbstractArray{<:Any,N}},
    others::SetCardProd{<:AbstractArray{<:Any,N}}...,
) where {N}
    SetCardProd(union.(s.sets, map(componentsets, others)...))
end

componentsets(s::SetCardProd) = s.sets

setcartprod(sets::AbstractArray{<:SetLike}) = SetCardProd(sets)
setcartprod(sets::Tuple{Vararg{SetLike}}) = SetCardProd(sets)
setcartprod(sets::NamedTuple{names,<:Tuple{Vararg{SetLike}}}) = SetCardProd(sets)

"""
    struct SetCartPower <: AbstractCartSetProd

Represents the n-fold Cartesian product of a set.
"""
struct SetCartPower{S,A} <: AbstractCartSetProd
    parent::S
    axes::A
end

maybestatic_length(s::SetCartPower) = size2length(maybestatic_size(s))
maybestatic_size(s::SetCartPower) = axes2size(s.axes)

@inline setcartpower(s::SetLike, dims) = SetCartPower(s, asaxes(dims))

componentsets(d::SetCartPower) = fill_with(d.parent, d.axes)

function Base.in(x::AbstractArray, s::SetCartPower)
    axes2size(s.axes) == size(x) ||
        throw(ArgumentError("Size of SetCartPower and given point are incompatible."))
    all(Base.Fix1(in, s.parent), x)
end

Base.isempty(s::SetCartPower) = isempty(s.parent) || size2length(axes2size(s.axes)) == 0

function Base.union(s::SetCartPower, others::SetCartPower...)
    axs = s.axes

    all(isequal(axs), map(x -> x.axes, others)) || throw(
        ArgumentError("Cannot create union of SetCartPower sets with different axes."),
    )

    setcartpower(union(s.parent, map(x -> x.parent, others)...))
end
