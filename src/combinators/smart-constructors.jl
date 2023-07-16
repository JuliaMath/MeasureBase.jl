
# Canonical measure type nesting, outer to inner:
#
# WeightedMeasure, Dirac, PowerMeasure, ProductMeasure


###############################################################################
# Half

half(μ::AbstractMeasure) = Half(μ)

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

@inline powermeasure(μ, exponent) = _generic_powermeasure_impl(asmeasure(μ), _pm_axes(exponent))

@inline _generic_powermeasure_stage1(μ::AbstractMeasure, ::Tuple{}) = μ

@inline function _generic_powermeasure_stage1(μ::AbstractMeasure, exponent::Tuple)
    _generic_powermeasure_stage2(μ, exponent)
end

@inline _generic_powermeasure_stage2(μ::AbstractMeasure, exponent::Tuple) = PowerMeasure(μ, exponent)

@inline function _generic_powermeasure_stage2(μ::Dirac, exponent::Tuple)
    Dirac(maybestatic_fill(μ.value, exponent))
end

@inline function _generic_powermeasure_stage2(μ::WeightedMeasure, exponent::Tuple)
    ν = μ.base^exponent
    k = maybestatic_length(ν) * μ.logweight
    return weightedmeasure(k, ν)
end


###############################################################################
# ProductMeasure

"""
    productmeasure(μs)

Constructs a product over a collection `μs` of measures.

Examples:

```julia
using MeasureBase, AffineMaps
productmeasure((StdNormal(), StdExponential()))
productmeasure(a = StdNormal(), b = StdExponential()))
productmeasure([pushfwd(Mul(scale), StdExponential()) for scale in 0.1:0.2:2])
productmeasure((pushfwd(Mul(scale), StdExponential()) for scale in 0.1:0.2:2))
"""
function productmeasure end
export productmeasure

@inline productmeasure(mar) = _generic_procuctmeasure_impl(mar)

@inline _generic_procuctmeasure_impl(mar::Fill) = powermeasure(_fill_value(mar), _fill_axes(mar))

@inline _generic_procuctmeasure_impl(mar::Tuple{Vararg{AbstractMeasure}}) = ProductMeasure(mar)
_generic_procuctmeasure_impl(mar::Tuple{Vararg{Dirac}}) = Dirac(map(m -> m.value), mar)
_generic_procuctmeasure_impl(mar::Tuple) = productmeasure(map(asmeasure, mar))

@inline _generic_procuctmeasure_impl(mar::NamedTuple{names,<:Tuple{Vararg{AbstractMeasure}}}) where names = ProductMeasure(mar)
_generic_procuctmeasure_impl(mar::NamedTuple{names,<:Tuple{Vararg{Dirac}}}) where names = Dirac(map(m -> m.value), mar)
_generic_procuctmeasure_impl(mar::NamedTuple) = productmeasure(map(asmeasure, mar))

@inline _generic_procuctmeasure_impl(mar::AbstractArray{<:AbstractProductMeasure}) = ProductMeasure(mar)
_generic_procuctmeasure_impl(mar::AbstractArray{<:Dirac}) = Dirac((m -> m.value).(mar))
_generic_procuctmeasure_impl(mar::AbstractArray) = ProductMeasure(asmeasure.(mar))

@inline function _generic_procuctmeasure_impl(mar::ReadonlyMappedArray{T,N,A,Returns{M}}) where {T,N,A,M}
    return powermeasure(mar.f.value, axes(mar.data))
end

@inline _generic_procuctmeasure_impl(mar::Base.Generator) = ProductMeasure(mar)

# TODO: Make this static when its length is static
@inline function _generic_procuctmeasure_impl(
    mar::AbstractArray{<:WeightedMeasure{StaticFloat64{W},M}},
) where {W,M}
    return weightedmeasure(W * length(mar), productmeasure(map(basemeasure, mar)))
end

# ToDo: Remove or at least refactor this (ProductMeasure shouldn't take a kernel at it's argument).

productmeasure(f, param_maps, pars) = ProductMeasure(kernel(f, param_maps), pars)

function productmeasure(k::ParameterizedTransitionKernel, pars)
    productmeasure(k.suff, k.param_maps, pars)
end

function productmeasure(f::Returns{W}, ::typeof(identity), pars) where {W<:WeightedMeasure}
    ℓ = _logweight(f.value)
    base = basemeasure(f.value)
    newbase = productmeasure(Returns(base), identity, pars)
    weightedmeasure(length(pars) * ℓ, newbase)
end

###############################################################################
# RestrictedMeasure
export restrict

restrict(f, b) = RestrictedMeasure(f, b)

###############################################################################
# SuperpositionMeasure

superpose(a::AbstractArray) = SuperpositionMeasure(a)

superpose(t::Tuple) = SuperpositionMeasure(t)
superpose(nt::NamedTuple) = SuperpositionMeasure(nt)

function superpose(μ::T, ν::T) where {T<:AbstractMeasure}
    if μ == ν
        return weightedmeasure(static(float(logtwo)), μ)
    else
        return superpose((μ, ν))
    end
end

function superpose(μ::AbstractMeasure, μs...)
    if all(==(μ), μs)
        return weightedmeasure(log(length(μs) + 1), μ)
    else
        return superpose((μ, μs...))
    end
end

add_measures(μs::AbstractVector, νs) = push!(μs, νs...)
add_measures(μs::Tuple, νs) = (μs..., νs...)

function superpose(μ::SuperpositionMeasure, μs...)
    SuperpositionMeasure(add_measures(μ.components, μs))
end

superpose(μ::SuperpositionMeasure) = μ

###############################################################################
# WeightedMeasure

function weightedmeasure(ℓ::R, b::M) where {R,M}
    WeightedMeasure{R,M}(ℓ, b)
end

function weightedmeasure(ℓ, b::WeightedMeasure)
    weightedmeasure(ℓ + _logweight(b), b.base)
end

###############################################################################
# TransitionKernel

# kernel(Normal(μ=2))
function kernel(μ::M) where {M<:ParameterizedMeasure}
    kernel(M)
end

function kernel(d::PowerMeasure)
    Base.Fix2(powermeasure, d.axes) ∘ kernel(d.parent)
end

function kernel(f)
    T = Core.Compiler.return_type(f, Tuple{Any})
    _kernel(f, T)
end

function _kernel(f, ::Type{T}) where {T}
    GenericTransitionKernel(f)
end

function _kernel(f, ::Type{P}) where {N,P<:ParameterizedMeasure{N}}
    k = length(N)
    C = constructorof(P)
    maps = ntuple(Val(k)) do i
        x -> @inbounds x[i]
    end

    kernel(params ∘ f, C, NamedTuple{N}(maps))
end

kernel(f::F, ::Type{M}; kwargs...) where {F<:Function,M} = kernel(f, M, NamedTuple(kwargs))

function kernel(f::F, ::Type{M}, nt::NamedTuple) where {F<:Function,M}
    ParameterizedTransitionKernel(M, f, nt)
end

function kernel(f::F, ::Type{M}, ::NamedTuple{()}) where {F<:Function,M}
    T = Core.Compiler.return_type(f, Tuple{Any})
    _kernel(f, M, T)
end

kernel(::Type{P}, nt::NamedTuple) where {P<:ParameterizedMeasure} = kernel(identity, P, nt)

kernel(::Type{T}; kwargs...) where {T} = kernel(T, NamedTuple(kwargs))

function kernel(::Type{M}, ::NamedTuple{()}) where {M}
    C = constructorof(M)
    TypedTransitionKernel(C, identity)
end

function _kernel(f::F, ::Type{M}, ::Type{NT}) where {M,F,N,NT<:NamedTuple{N}}
    k = length(N)
    maps = ntuple(Val(k)) do i
        x -> @inbounds x[i]
    end

    ParameterizedTransitionKernel(M, values ∘ f, NamedTuple{N}(maps))
end

kernel(f::F; kwargs...) where {F<:Function} = kernel(f, NamedTuple(kwargs))

function kernel(f::F, nt::NamedTuple{()}) where {F<:Function}
    T = Core.Compiler.return_type(f, Tuple{Any})
    _kernel(f, T)
end
