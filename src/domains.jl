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
    const MeasureBase.SetLike = Union{MeasureBase.ValueSet, Base.AbstractSet, IntervalSets.Domain}

Any kind of (measurable) set.

There needs to be an implicit sigma-algebra for subtypes of
`MeasureBase.SetLike` to make them useable for measures. This can't easily be
imposed via type constraints, though, so it is by-contract.
"""
const SetLike = Union{MeasureBase.ValueSet,Base.AbstractSet,IntervalSets.Domain}

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
    RealValues() isa MeasureBase.ValueSet

The real numbers.
"""
struct RealValues <: ValueSet end

@inline Base.in(x::Real, ::RealValues) = true
@inline Base.in(x, ::RealValues) = isreal(x)

@inline Base.isempty(::RealValues) = false

@inline Base.union(s::RealValues, ::RealValues...) = s

@inline Base.minimum(::RealValues) = static(-Inf)
@inline Base.maximum(::RealValues) = static(Inf)

testvalue(::Type{T}, ::RealValues) where {T} = zero(T)

"""
    const MeasureBase.ℝ = RealValues()

The set of all real numbers, see [`MeasureBase.RealValues`](@ref).
"""
const ℝ = RealValues()
export ℝ

Base.show(io::IO, ::RealValues) = print(io, "ℝ")
Base.show(io::IO, ::MIME"text/plain", ::RealValues) = print(io, "MeasureBase.ℝ")

"""
    MeasureBase.IntegerValues() isa MeasureBase.ValueSet
"""
struct IntegerValues <: ValueSet end

@inline Base.in(x::Integer, ::IntegerValues) = true
@inline Base.in(x, ::IntegerValues) = isinteger(x)

@inline Base.isempty(::IntegerValues) = false

@inline Base.union(s::IntegerValues, ::IntegerValues...) = s

testvalue(::Type{T}, ::IntegerValues) where {T} = zero(T)

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
export ℤ

Base.show(io::IO, ::IntegerValues) = print(io, "ℤ")
Base.show(io::IO, ::MIME"text/plain", ::IntegerValues) = print(io, "MeasureBase.ℤ")

"""
    struct MeasureBase.BoundedInts{L,U} <: MeasureBase.ValueSet

The integers from `lower` to `upper` (bounds may be infinite).

Constructors:

```julia
BoundedInts(lower, upper)
ℤ[lower:upper]
```
"""
struct BoundedInts{L,U} <: ValueSet
    lower::L
    upper::U
end

@inline Base.in(x, b::BoundedInts) = x ∈ ℤ && b.lower <= x <= b.upper

Base.isempty(b::BoundedInts) = b.lower > b.upper

Base.minimum(b::BoundedInts) = b.lower
Base.maximum(b::BoundedInts) = b.upper

function Base.show(io::IO, b::BoundedInts)
    io = IOContext(io, :compact => true)
    print(io, "ℤ[", b.lower, ":", b.upper, "]")
end

testvalue(b::BoundedInts) = convert(Int, clamp(0, dynamic(b.lower), dynamic(b.upper)))
testvalue(::Type{T}, b::BoundedInts) where {T} = convert(T, testvalue(b))

Base.getindex(::typeof(ℤ), r::AbstractUnitRange) = BoundedInts(extrema(r)...)

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
setcartprod(sets::NamedTuple{names,<:Tuple{Vararg{SetLike}}}) where {names} =
    CartesianProduct(sets)

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

@inline Base.isempty(s::CartesianProduct) = any(isempty, componentsets(s))

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
@inline pwr_size(s::CartesianPower) = axes2size(pwr_axes(s))

componentsets(s::CartesianPower) = maybestatic_fill(pwr_base(s), pwr_axes(s))

function Base.in(x::AbstractArray, s::CartesianPower)
    pwr_size(s) == size(x) ||
        throw(ArgumentError("Size of CartesianPower and given point are incompatible."))
    isempty(x) ? true : all(Base.Fix2(in, pwr_base(s)), x)::Bool
end

Base.isempty(s::CartesianPower) = isempty(pwr_base(s)) || size2length(pwr_size(s)) == 0

function Base.union(s::CartesianPower, others::CartesianPower...)
    axs = pwr_axes(s)

    all(isequal(axs), map(pwr_axes, others)) || throw(
        ArgumentError("Cannot create union of CartesianPower sets with different axes."),
    )

    setcartpower(union(pwr_base(s), map(pwr_base, others)...), axs)
end


"""
    struct CombinedSet <: ValueSet

Represents a combination of two sets.

User code should not create instances of `CombinedSet` directly, but should
call [`combinesets(f_c, α, β)`](@ref) instead.
"""
struct CombinedSet{FC,MA<:SetLike,MB<:SetLike} <: ValueSet
    f_c::FC
    α::MA
    β::MB
end

function Base.in(@nospecialize(x), ::CombinedSet)
    throw(ArgumentError("Cannot test if a value lies within a combined set."))
end

maybe_in(@nospecialize(x), ::CombinedSet) = true

Base.isempty(s::CombinedSet) = isempty(s.α) || isempty(s.β)

"""
    combinesets(f_c, α, β)

Combine two sets `α` and `β` into the set of all values `f_c(a, b)` with
`a ∈ α` and `b ∈ β`.

`f_c` must combine values as described in [`mcombine`](@ref). Uses set
representations more specific than [`MeasureBase.CombinedSet`](@ref) where
possible.
"""
function combinesets end
export combinesets

@inline combinesets(f_c, α::SetLike, β::SetLike) = _generic_combinesets(f_c, α, β)

# Combining the implicit domains of two measures yields the implicit domain
# of the combined measure:
@inline combinesets(f_c, α::ImplicitDomain, β::ImplicitDomain) =
    ImplicitDomain(mcombine(f_c, α.m, β.m))

@inline _generic_combinesets(::typeof(firstarg), α::SetLike, β::SetLike) = α
@inline _generic_combinesets(::typeof(secondarg), α::SetLike, β::SetLike) = β
@inline _generic_combinesets(::typeof(tuple), α::SetLike, β::SetLike) =
    setcartprod((α, β))
@inline _generic_combinesets(f_c::typeof(vcat), α::SetLike, β::SetLike) =
    _combinesets_cat(f_c, α, β)
@inline _generic_combinesets(f_c::typeof(merge), α::SetLike, β::SetLike) =
    _combinesets_cat(f_c, α, β)
@inline _generic_combinesets(f_c, α::SetLike, β::SetLike) = CombinedSet(f_c, α, β)

_combinesets_cat(
    ::typeof(vcat),
    α::CartesianProduct{<:AbstractVector},
    β::CartesianProduct{<:AbstractVector},
) = setcartprod(vcat(componentsets(α), componentsets(β)))

_combinesets_cat(
    ::typeof(merge),
    α::CartesianProduct{<:NamedTuple},
    β::CartesianProduct{<:NamedTuple},
) = setcartprod(merge(componentsets(α), componentsets(β)))

# Concatenating one-dimensional powers of equal base sets yields a longer
# power. Set equality can typically only be established at runtime, so this
# simplification only happens when it is decidable from the set types alone:
function _combinesets_cat(
    ::typeof(vcat),
    α::CartesianPower{<:Any,<:Tuple{Any}},
    β::CartesianPower{<:Any,<:Tuple{Any}},
)
    if _static_isequal(pwr_base(α), pwr_base(β)) isa True
        n = size2length(pwr_size(α)) + size2length(pwr_size(β))
        setcartpower(pwr_base(α), (n,))
    else
        CombinedSet(vcat, α, β)
    end
end

_combinesets_cat(f_c, α::SetLike, β::SetLike) = CombinedSet(f_c, α, β)
