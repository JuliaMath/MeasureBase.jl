# Custom abstract set type. Design reserve to be able to switch to
#`Base.AbstractSet` or another set type hierarchy in the future:
"""
    MeasureBase.ValueSet

Abstract type for some measurable sets.

Not every measurable set needs to be a a subtype of
`MeasureBase.ValueSet`.

See also [`MeasureBase.SetLike`](@ref).
"""
abstract type ValueSet end

"""
    const MeasureBase.SetLike = Union{MeaureBase.ValueSet, Base.AbstractSet, IntervalSets.Domain}

Any kind of (measurable) set.

There needs to be an implicit sigma-algebra for subtypes of
`MeasureBase.SetLike` to make them useable for measures. This can't easily be
imposed via type constraints, though, to is is by-contract.
"""
const SetLike = Union{MeasureBase.ValueSet,Base.AbstractSet,IntervalSets.Domain}

"""
    mdomain(m)::MeasureBase.SetLike

Return the domain, i.e. the measurable set, of the measure `m`.

The measure must allow for evaluating densities and the like over the whole
domain, even if the support of the measure is only a subset of the domain.

May return [`MeasureBase.ImplicitDomain`](@ref) if the domain cannot be
computed (efficiently).
"""
function mdomain end
export mdomain

@inline mdomain(m) = ImplicitDomain(m)

"""
    valdomain(x)::MeasureBase.SetLike

Return the domain of a given value.

May return [`MeasureBase.UnknownDomain`](@ref) if no domain type is available
that can represents values like `x`.
"""
function valdomain end
export valdomain

@inline valdomain(x) = UnknownDomain(x)

"""
    MeasureBase.maybe_in(x, s)

Test if `x` may be a member of `s`.

Defaults to `in(x, s)`, but may be specialized for certain types of `s`,
e.g. for `s::MeasureBase.ImplicitDomain`.
"""
function maybe_in end

maybe_in(x, s) = in(x, s)

"""
    struct MeasureBase.ImplicitDomain{M} <: MeasureBase.ValueSet

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
struct ImplicitDomain{M} <: ValueSet
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
    struct MeasureBase.UnknownDomain{T} <: MeasureBase.ValueSet

Represents the unknown domain of a value of type `T`.

Constructors:

```
MeasureBase.UnknownDomain(x::T)
```

Does not support `Base.in(x, s::MeasureBase.UnknownDomain)`, and
`MeasureBase.maybe_in(x, s::MeasureBase.UnknownDomain)` always return `true`
(unless specialized for the measure type).

`isempty` will always return false, `UnknownDomain` should only be created
if a value of type `T` existed in the first place, which implies that the
domain can not be empty.
"""
struct UnknownDomain{T} <: ValueSet end

UnknownDomain(::T) where {T} = UnknownDomain{T}()

Base.eltype(::UnknownDomain{T}) where {T} = T

@inline Base.union(s::UnknownDomain, others::UnknownDomain...) =
    UnknownDomain{promote_type(eltype(s), map(eltype, others)...)}()

function Base.in(@nospecialize(x), ::UnknownDomain)
    throw(ArgumentError("Cannot test if a value lies withing an unknown domain."))
end

maybe_in(@nospecialize(x), ::UnknownDomain) = true

Base.isempty(::UnknownDomain) = false

"""
    RealInterval() isa MeasureBase.ValueSet

The real numbers.
"""
struct RealValues <: ValueSet end

@inline Base.in(x::Real, ::RealValues) = true
@inline Base.in(x, ::RealValues) = isreal(x)

@inline Base.isempty(::RealValues) = false

@inline Base.union(s::RealValues, ::RealValues...) = s

@inline Base.minimum(::RealValues) = static(-Inf)
@inline Base.maximum(::RealValues) = static(Inf)

"""
    const MeasureBase.ℝ = RealValues()

The set of all real numbers, see [`MeasureBase.RealValues`](@ref).
"""
const ℝ = RealValues()

Base.show(io::IO, ::MIME"text/plain", ::RealValues) = print(io, "MeasureBase.ℝ")

"""
    MeasureBase.IntegerValues() isa MeasureBase.ValueSet
"""
struct IntegerValues <: ValueSet end

@inline Base.in(x::Integer, ::IntegerValues) = true
@inline Base.in(x, ::IntegerValues) = isinteger(x)

@inline Base.isempty(::IntegerValues) = false

@inline Base.union(s::IntegerValues, ::IntegerValues...) = s

# # This could get tricky with mixed-precision code. Probably needs some
# # special AbstractInteger infinity type (but custom AbstractInteger types
# # may cause a lot of method invalidations, which is why Static.StaticInteger
# # is not an AbstractInteger).
# @inline Base.minimum(::RealValues) = static(typemax(Int64))
# @inline Base.maximum(::RealValues) = static(typemin(Int64))

"""
    const ℤ = IntegerValues()

The set of all integers, see [`MeasureBase.IntegerValues`](@ref).
"""
const ℤ = IntegerValues()

Base.show(io::IO, ::MIME"text/plain", ::IntegerValues) = print(io, "MeasureBase.ℤ")

"""
    struct MeasureBase.AbstractCartSetProd <: ValueSet

Supertype for cartesian products of sets.
"""
abstract type AbstractCartSetProd <: ValueSet end

"""
    struct CartesianProduct <: AbstractCartSetProd

A cartesian product over a collection of sets.

Constructor:

```julia
prodset = CartesianProduct(sets)
```

`sets` may be a `Tuple`, `NamedTuple` or `AbstractArray` of sets/domains.
"""
struct CartesianProduct{S<:Union{Tuple,NamedTuple,AbstractArray}} <: AbstractCartSetProd
    _sets::S
end

componentsets(s::CartesianProduct) = s._sets

setcartprod(sets::AbstractArray{<:SetLike}) = CartesianProduct(sets)
setcartprod(sets::Tuple{Vararg{SetLike}}) = CartesianProduct(sets)
setcartprod(sets::NamedTuple{names,<:Tuple{Vararg{SetLike}}}) = CartesianProduct(sets)

@inline Base.in(x::Tuple{}, s::CartesianProduct{Tuple{}}) = true
@inline Base.in(x::Tuple{Vararg{Any,N}}, s::CartesianProduct{<:Tuple{Vararg{Any,N}}}) where {N} =
    prod(map(in, x, componentsets(s)))::Bool
@inline Base.in(x::NamedTuple{names}, s::CartesianProduct{<:NamedTuple{names}}) where {names} =
    prod(map(in, values(x), values(componentsets(s))))::Bool
# ToDo: Allow this?
# Base.in(x::AbstractVector, s::CartesianProduct{<:Tuple}) = all(in.(x,componentsets(s)))::Bool
function Base.in(
    x::AbstractArray{<:Any,N},
    s::CartesianProduct{<:AbstractArray{<:Any,N}},
) where {N}
    sets = componentsets(s)
    isempty(x) && isempty(sets) ? true : all(in.(x, sets))::Bool
end

@inline Base.isempty(s::CartesianProduct) = all(!isempty, componentsets(s))

@inline function Base.union(
    s::CartesianProduct{<:Tuple{Vararg{Any,N}}},
    others::CartesianProduct{<:Tuple{Vararg{Any,N}}}...,
) where {N}
    CartesianProduct(map(union, componentsets(s), map(componentsets, others)...))
end

@inline function Base.union(
    s::CartesianProduct{<:NamedTuple{names}},
    others::CartesianProduct{<:NamedTuple{names}}...,
) where {names}
    CartesianProduct(map(union, componentsets(s), map(componentsets, others)...))
end

function Base.union(
    s::CartesianProduct{<:AbstractArray{<:Any,N}},
    others::CartesianProduct{<:AbstractArray{<:Any,N}}...,
) where {N}
    CartesianProduct(union.(componentsets(s), map(componentsets, others)...))
end

"""
    struct CartesianPower <: AbstractCartSetProd

Represents the n-fold Cartesian product of a set.
"""
struct CartesianPower{S,A} <: AbstractCartSetProd
    _base::S
    _axes::A
end

@inline setcartpower(s::SetLike, dims) = CartesianPower(s, asaxes(dims))

@inline pwr_base(s::CartesianPower) = s._base
@inline pwr_axes(s::CartesianPower) = s._axes
@inline pwr_size(s::CartesianPower) = axes2size(s.axes)

componentsets(d::CartesianPower) = fill_with(d.parent, d.axes)

function Base.in(x::AbstractArray, s::CartesianPower)
    axes2size(s.axes) == size(x) ||
        throw(ArgumentError("Size of CartesianPower and given point are incompatible."))
    isempty(x) ? true : all(Base.Fix1(in, s.parent), x)::Bool
end

Base.isempty(s::CartesianPower) = isempty(s.parent) || size2length(axes2size(s.axes)) == 0

function Base.union(s::CartesianPower, others::CartesianPower...)
    axs = s.axes

    all(isequal(axs), map(x -> x.axes, others)) || throw(
        ArgumentError("Cannot create union of CartesianPower sets with different axes."),
    )

    setcartpower(union(s.parent, map(x -> x.parent, others)...))
end
