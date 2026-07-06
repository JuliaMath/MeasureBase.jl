
# Canonical measure type nesting, outer to inner:
#
# WeightedMeasure, Dirac, PowerMeasure, ProductMeasure


###############################################################################
# Half

"""
    half(μ::AbstractMeasure)

Constructs the half-measure of a measure `μ` that is symmetric around zero:
`μ` folded onto the non-negative half-line.
"""
half(μ::AbstractMeasure) = Half(μ)
export half

###############################################################################
# PowerMeaure

"""
    powermeasure(μ, dims)
    powermeasure(μ, axes)

Constructs a power of a measure `μ`.

`powermeasure(μ, exponent)` is semantically equivalent to
`productmeasure(Fill(μ, exponent))`, but more efficient.
"""
function powermeasure end
export powermeasure

@inline powermeasure(μ, exponent) = _generic_powermeasure_stage1(asmeasure(μ), asaxes(exponent))

@inline _generic_powermeasure_stage1(μ::AbstractMeasure, ::Tuple{}) = μ

@inline function _generic_powermeasure_stage1(μ::AbstractMeasure, exponent::Tuple)
    _generic_powermeasure_stage2(μ, exponent)
end

@inline _generic_powermeasure_stage2(μ::AbstractMeasure, exponent::Tuple) =
    PowerMeasure(μ, exponent)

@inline function _generic_powermeasure_stage2(μ::Dirac, exponent::Tuple)
    Dirac(maybestatic_fill(μ.x, exponent))
end

@inline function _generic_powermeasure_stage2(μ::WeightedMeasure, exponent::Tuple)
    ν = μ.base^exponent
    k = size2length(axes2size(exponent)) * μ.logweight
    return weightedmeasure(k, ν)
end

###############################################################################
# ProductMeasure

"""
    productmeasure(μs)

Constructs a product over a collection `μs` of measures.

Examples:

```julia
productmeasure((StdNormal(), StdExponential()))
productmeasure((a = StdNormal(), b = StdExponential()))
productmeasure([pushfwd(Base.Fix1(*, scale), StdExponential()) for scale in 0.1:0.2:2])
```
"""
function productmeasure end
export productmeasure

@inline productmeasure(mar) = _generic_productmeasure_impl(mar)

@inline _generic_productmeasure_impl(mar::FillArrays.Fill) =
    powermeasure(_fill_value(mar), _fill_axes(mar))

# Empty products are unit measures:
@inline _generic_productmeasure_impl(::Tuple{}) = Dirac(())
@inline _generic_productmeasure_impl(::NamedTuple{()}) = Dirac(NamedTuple())

@inline _generic_productmeasure_impl(mar::Tuple{Vararg{AbstractMeasure}}) =
    ProductMeasure(mar)
_generic_productmeasure_impl(mar::Tuple{Dirac,Vararg{Dirac}}) = Dirac(map(m -> m.x, mar))
_generic_productmeasure_impl(mar::Tuple{WeightedMeasure,Vararg{WeightedMeasure}}) =
    weightedmeasure(sum(map(_logweight, mar)), productmeasure(map(m -> m.base, mar)))
_generic_productmeasure_impl(mar::Tuple) = productmeasure(map(asmeasure, mar))

@inline _generic_productmeasure_impl(
    mar::NamedTuple{names,<:Tuple{Vararg{AbstractMeasure}}},
) where {names} = ProductMeasure(mar)
_generic_productmeasure_impl(
    mar::NamedTuple{names,<:Tuple{Dirac,Vararg{Dirac}}},
) where {names} = Dirac(map(m -> m.x, mar))
_generic_productmeasure_impl(
    mar::NamedTuple{names,<:Tuple{WeightedMeasure,Vararg{WeightedMeasure}}},
) where {names} =
    weightedmeasure(sum(map(_logweight, values(mar))), productmeasure(map(m -> m.base, mar)))
_generic_productmeasure_impl(mar::NamedTuple) = productmeasure(map(asmeasure, mar))

_generic_productmeasure_impl(mar::AbstractArray{<:Dirac}) = Dirac((m -> m.x).(mar))

_generic_productmeasure_impl(mar::AbstractArray{<:WeightedMeasure}) =
    weightedmeasure(sum(_logweight, mar), productmeasure((m -> m.base).(mar)))

@inline function _generic_productmeasure_impl(
    mar::AbstractArray{<:WeightedMeasure{StaticFloat64{W},M}},
) where {W,M}
    return weightedmeasure(
        static(W) * maybestatic_length(mar),
        productmeasure((m -> m.base).(mar)),
    )
end

function _generic_productmeasure_impl(mar::AbstractArray{T}) where {T}
    if Base.issingletontype(T)
        powermeasure(instance(T), axes(mar))
    elseif T <: AbstractMeasure
        ProductMeasure(mar)
    else
        ProductMeasure(map(asmeasure, mar))
    end
end

@inline function _generic_productmeasure_impl(
    mar::ReadonlyMappedArray{T,N,A,Returns{M}},
) where {T,N,A,M}
    return powermeasure(mar.f.value, axes(mar.data))
end

@inline _generic_productmeasure_impl(mar::Base.Generator) = ProductMeasure(mar)

###############################################################################
# PushforwardMeasure

# The pushforward of a point mass is a point mass. Note that no density
# volume correction applies, Dirac measures are density-defined relative
# to counting measure:
_pushfwd_impl(f, μ::Dirac, ::PushFwdStyle) = Dirac(f(μ.x))

# Pushforward and weighting commute:
function _pushfwd_impl(f, μ::WeightedMeasure, style::PushFwdStyle)
    weightedmeasure(μ.logweight, _pushfwd_impl(f, μ.base, style))
end

###############################################################################
# RestrictedMeasure
export restrict

@inline restrict(f) = Base.Fix1(restrict, f)

restrict(f, μ) = RestrictedMeasure(f, asmeasure(μ))

# Nested restrictions fuse into a single predicate:
restrict(f, μ::RestrictedMeasure) =
    RestrictedMeasure(x -> μ.predicate(x) && f(x), μ.base)

###############################################################################
# SuperpositionMeasure

"""
    superpose(μs...)
    superpose(μs)

Constructs a superposition of measures, given either as separate arguments or
as a collection (array, tuple or named tuple) of measures.

The vararg form simplifies algebraically: equal measures combine into weighted
measures (`superpose(μ, μ) == weightedmeasure(log(2), μ)`), weighted measures
with equal bases add their weights, and superpositions merge their components.
To keep `superpose` type stable, such simplifications only happen when
equality of the measures involved can be decided from their types alone.
Collections are wrapped as-is, apart from cost-free structural simplifications.
"""
function superpose end
export superpose

superpose(μ::AbstractMeasure) = μ

function superpose(μ::AbstractMeasure, ν::AbstractMeasure, more::AbstractMeasure...)
    superpose(_superpose_two(μ, ν), more...)
end

function superpose(a::AbstractArray{T}) where {T}
    if Base.issingletontype(T)
        weightedmeasure(log(length(a)), asmeasure(instance(T)))
    else
        SuperpositionMeasure(a)
    end
end

superpose(a::FillArrays.Fill) = weightedmeasure(log(length(a)), asmeasure(_fill_value(a)))

superpose(t::Tuple) = SuperpositionMeasure(t)
superpose(nt::NamedTuple) = SuperpositionMeasure(nt)

# Measure equality can typically only be established at runtime, but measure
# construction must be type stable, so simplifications may only depend on
# measure equality that is decidable from the measure types alone:
@inline _static_isequal(::T, ::T) where {T} = static(Base.issingletontype(T))
@inline _static_isequal(::Any, ::Any) = static(false)

function _superpose_two(μ::AbstractMeasure, ν::AbstractMeasure)
    if _static_isequal(μ, ν) isa True
        weightedmeasure(static(float(logtwo)), μ)
    else
        SuperpositionMeasure((μ, ν))
    end
end

function _superpose_two(μ::WeightedMeasure, ν::WeightedMeasure)
    if _static_isequal(μ.base, ν.base) isa True
        weightedmeasure(logaddexp(asnonstatic(μ.logweight), asnonstatic(ν.logweight)), μ.base)
    else
        SuperpositionMeasure((μ, ν))
    end
end

function _superpose_two(μ::WeightedMeasure, ν::AbstractMeasure)
    if _static_isequal(μ.base, ν) isa True
        weightedmeasure(log1pexp(asnonstatic(μ.logweight)), μ.base)
    else
        SuperpositionMeasure((μ, ν))
    end
end

function _superpose_two(μ::AbstractMeasure, ν::WeightedMeasure)
    if _static_isequal(μ, ν.base) isa True
        weightedmeasure(log1pexp(asnonstatic(ν.logweight)), ν.base)
    else
        SuperpositionMeasure((μ, ν))
    end
end

_superpose_two(μ::SuperpositionMeasure, ν::SuperpositionMeasure) =
    SuperpositionMeasure(_cat_measures(μ.components, ν.components))
_superpose_two(μ::SuperpositionMeasure, ν::AbstractMeasure) =
    SuperpositionMeasure(_cat_measures(μ.components, (ν,)))
_superpose_two(μ::AbstractMeasure, ν::SuperpositionMeasure) =
    SuperpositionMeasure(_cat_measures((μ,), ν.components))
_superpose_two(μ::SuperpositionMeasure, ν::WeightedMeasure) =
    SuperpositionMeasure(_cat_measures(μ.components, (ν,)))
_superpose_two(μ::WeightedMeasure, ν::SuperpositionMeasure) =
    SuperpositionMeasure(_cat_measures((μ,), ν.components))

###############################################################################
# WeightedMeasure

"""
    weightedmeasure(logweight::Real, μ)

Constructs a measure that behaves like the measure `μ`, but with its density
scaled by `exp(logweight)`. Weights of nested weighted measures combine
additively.
"""
function weightedmeasure end
export weightedmeasure

function weightedmeasure(ℓ::R, b::M) where {R,M}
    WeightedMeasure{R,M}(ℓ, b)
end

function weightedmeasure(ℓ, b::WeightedMeasure)
    weightedmeasure(ℓ + _logweight(b), b.base)
end

