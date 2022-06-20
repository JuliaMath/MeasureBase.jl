
###############################################################################
# Half

half(μ::AbstractMeasure) = Half(μ)

###############################################################################
# PointwiseProductMeasure

function pointwiseproduct(μ::AbstractMeasure, ℓ::Likelihood)
    T = Core.Compiler.return_type(ℓ.k, Tuple{gentype(μ)})
    return pointwiseproduct(T, μ, ℓ)
end

function pointwiseproduct(::Type{T}, μ::AbstractMeasure, ℓ::Likelihood) where {T}
    return PointwiseProductMeasure(μ, ℓ)
end

###############################################################################
# PowerMeaure

powermeasure(m::AbstractMeasure, ::Tuple{}) = m

function powermeasure(
    μ::WeightedMeasure,
    dims::Tuple{<:AbstractArray, Vararg{<:AbstractArray}},
)
    k = mapreduce(length, *, dims) * μ.logweight
    return weightedmeasure(k, μ.base^dims)
end

function powermeasure(μ::WeightedMeasure, dims::Tuple{<:Integer, Vararg{<:Integer}})
    k = prod(dims) * μ.logweight
    return weightedmeasure(k, μ.base^dims)
end

###############################################################################
# ProductMeasure

productmeasure(mar::Fill) = powermeasure(mar.value, mar.axes)

function productmeasure(mar::ReadonlyMappedArray{T,N,A,Returns{M}}) where {T,N,A,M}
    return powermeasure(mar.f.value, axes(mar.data))
end

productmeasure(mar::Base.Generator) = ProductMeasure(mar)
productmeasure(mar::AbstractArray) = ProductMeasure(mar)

# TODO: Make this static when its length is static
@inline function productmeasure(
    mar::AbstractArray{WeightedMeasure{StaticFloat64{W},M}},
) where {W,M}
    return weightedmeasure(W * length(mar), productmeasure(map(basemeasure, mar)))
end

productmeasure(nt::NamedTuple) = ProductMeasure(nt)
productmeasure(tup::Tuple) = ProductMeasure(tup)

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
        return weightedmeasure(logtwo, μ)
    else
        return superpose((μ, ν))
    end
end

function superpose(::T, ::T) where {T<:SuperpositionMeasure}
    @error "FIXME"
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

function kernel(::Type{M}, ::NamedTuple{()}) where {M<:ParameterizedMeasure}
    @error "FIXME"
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
