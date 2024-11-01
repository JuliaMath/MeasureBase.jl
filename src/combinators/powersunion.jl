
"""
    PowersUnionMeasure

A `ν = PowersUnionMeasure(μ)` represents the disjoint union, resp. the direct
sum, over all `μ^n` for all finite integers `n >= 0`.

`logdensity_of(ν, x)` requires x to be a vector or set of variates of `μ`.
Vectors are treated as unordered, so a vector `x` is implicitly treated as
`Set(x)`, and both result in the same log-density value.

`PowersUnionMeasure`s are primarily indended to be used as reference resp.
root measures for measures that represent random point processes and similar.
"""
struct PowersUnionMeasure{M} <: AbstractMeasure
    parent::M
end
export PowersUnionMeasure

const _SomeVectorOrSet{T} = Union{AbstractVector{T},AbstractSet{T}}

@inline basemeasure(pm::PowersUnionMeasure) = PowersUnionMeasure(basemeasure(pm.parent))
@inline rootmeasure(pm::PowersUnionMeasure) = PowersUnionMeasure(rootmeasure(pm.parent))

logdensity_def(pm::PowersUnionMeasure{M}, x) where {M} = sum(Base.Fix1(logdensity_def, pm.parent), x)
@inline logdensity_def(::PowersUnionMeasure{<:PrimitiveMeasure}, ::Any) = static(0.0)
@inline logdensity_def(::PowersUnionMeasure{<:PowerMeasure{<:PrimitiveMeasure}}, ::Any) = static(0.0)

# ToDo: Optimize for PowersUnionMeasure of primitive measures?
insupport(pm::PowersUnionMeasure, x::_SomeVectorOrSet) = all(dynamic ∘ Base.Fix1(insupport, pm.parent), x)
insupport(pm::PowersUnionMeasure, ::Any) = false

@inline getdof(pm::PowersUnionMeasure) = NoDOF(typeof(pm))

# ToDo: use checked_arg of parent measure if not a PowersUnionMeasure over a
# primitive measure? PowerMeasure doesn't check, currently.
@inline checked_arg(::PowersUnionMeasure, x::_SomeVectorOrSet) = x
checked_arg(::PowersUnionMeasure, ::Any) = throw(ArgumentError("Variates of PowersUnionMeasure must be vectors or sets"))

Base.rand(::AbstractRNG, ::Type, ::PowersUnionMeasure) = throw(ArgumentError("rand not available for PowersUnionMeasure"))
