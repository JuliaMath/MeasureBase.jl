
using LogarithmicNumbers
using LogExpFunctions

export SuperpositionMeasure

abstract type AbstractSuperpositionMeasure <: AbstractMeasure end

@doc raw"""
    struct SuperpositionMeasure{NT} <: AbstractMeasure
        components :: NT
    end
Superposition of measures is analogous to mixture distributions, but (because
measures need not be normalized) requires no scaling.
The superposition of two measures μ and ν can be more concisely written as μ + ν.
Superposition measures satisfy
    
    basemeasure(μ + ν) == basemeasure(μ) + basemeasure(ν)


```math
    \begin{aligned}\frac{\mathrm{d}(\mu+\nu)}{\mathrm{d}(\alpha+\beta)} & =\frac{f\,\mathrm{d}\alpha+g\,\mathrm{d}\beta}{\mathrm{d}\alpha+\mathrm{d}\beta}\\
     & =\frac{f\,\mathrm{d}\alpha}{\mathrm{d}\alpha+\mathrm{d}\beta}+\frac{g\,\mathrm{d}\beta}{\mathrm{d}\alpha+\mathrm{d}\beta}\\
     & =\frac{f}{1+\frac{\mathrm{d}\beta}{\mathrm{d}\alpha}}+\frac{g}{\frac{\mathrm{d}\alpha}{\mathrm{d}\beta}+1}\\
     & =\frac{f}{1+\left(\frac{\mathrm{d}\alpha}{\mathrm{d}\beta}\right)^{-1}}+\frac{g}{\frac{\mathrm{d}\alpha}{\mathrm{d}\beta}+1}\ .
    \end{aligned}
```
"""
struct SuperpositionMeasure{C} <: AbstractSuperpositionMeasure
    components::C
end

massof(m::SuperpositionMeasure) = sum(massof, m.components)

function Pretty.tile(d::SuperpositionMeasure)
    result = Pretty.literal("SuperpositionMeasure(")
    result *= Pretty.list_layout([Pretty.tile.(d.components)...])
    result *= Pretty.literal(")")
end

testvalue(::Type{T}, μ::SuperpositionMeasure) where {T} = testvalue(T, first(μ.components))

function Base.:+(μ::AbstractMeasure, ν::AbstractMeasure)
    superpose(μ, ν)
end

# Masks components outside of their support with -Inf:
@inline _masked_logd(ℓ, ins) = ifelse(_insupport_mask(ins), ℓ, oftype(ℓ, -Inf))

# Branch-free logsumexp over the components, valid for infinite entries:
@inline function _logsumexp_components(ℓs)
    m = reduce(max, ℓs)
    m_finite = ifelse(isfinite(m), m, zero(m))
    m_finite + log(sum(map(ℓ -> exp(ℓ - m_finite), ℓs)))
end

# The density of a superposition relative to the superposition of the
# component base measures, in log space: each component contributes its
# own density, divided by the density of the superposed base measures
# relative to its own base measure.
function logdensity_def(s::SuperpositionMeasure, x)
    cs = values(s.components)
    αs = map(basemeasure, cs)
    terms = map(cs, αs) do cᵢ, αᵢ
        ℓᵢ = _dynamic_logd(logdensity_def(cᵢ, x), x)
        log_dΣα_dαᵢ = _logsumexp_components(map(αs) do αⱼ
            _masked_logd(logdensity_rel(αⱼ, αᵢ, x), insupportof(αⱼ, x))
        end)
        _masked_logd(ℓᵢ - log_dΣα_dαᵢ, insupportof(cᵢ, x))
    end
    _logsumexp_components(terms)
end

function _superpos_logdensity_rel(s::SuperpositionMeasure, β, x)
    cs = values(s.components)
    ds = map(cs) do μ
        _masked_logd(logdensity_rel(μ, β, x), insupportof(μ, x))
    end
    _logsumexp_components(ds)
end

@inline logdensity_rel_def(s::SuperpositionMeasure, β, x) = _superpos_logdensity_rel(s, β, x)

@inline logdensity_rel_def(s::SuperpositionMeasure, β::SuperpositionMeasure, x) =
    _superpos_logdensity_rel(s, β, x)

@inline logdensity_rel_def(s, β::SuperpositionMeasure, x) = -_superpos_logdensity_rel(β, s, x)

@inline density_def(s::SuperpositionMeasure, x) = exp(logdensity_def(s, x))

# Pushforward is linear, superpositions push forward component by component,
# so that each component keeps or adapts its own root measure. Pushforwards
# by transports need the whole origin (see `basemeasure` of pushforwards):
function _pushfwd_impl(f, μ::SuperpositionMeasure, style::PushFwdStyle, varsize)
    superpose(map(c -> _pushfwd_impl(f, c, style, varsize), μ.components))
end
function _pushfwd_impl(f::TransportFunction, μ::SuperpositionMeasure, style::PushFwdStyle, varsize)
    _pushfwd_with_style(f, μ, _pushfwd_style(style, f, μ), varsize)
end

function basemeasure(μ::SuperpositionMeasure{<:Tuple})
    superpose(map(basemeasure, μ.components)...)
end

function basemeasure(μ::SuperpositionMeasure{<:AbstractArray})
    bases = map(basemeasure, μ.components)
    allequal(bases) ? weightedmeasure(log(length(bases)), first(bases)) : superpose(bases)
end

basemeasure(μ::SuperpositionMeasure) = superpose(map(basemeasure, μ.components))

# Superpositions of the local measures of the components (e.g. of binds):
localmeasure(μ::SuperpositionMeasure{<:Union{Tuple,NamedTuple}}, x) = SuperpositionMeasure(map(Base.Fix2(localmeasure, x), μ.components))

function _component_masses(μ::SuperpositionMeasure)
    masses = map(massof, values(μ.components))
    total = sum(masses)
    total isa AbstractUnknownMass && throw(
        ArgumentError("Cannot sample from a superposition of measures of unknown mass"),
    )
    return map(dynamic, masses), dynamic(total)
end

# Single variates of superpositions with a variate rank are drawn like
# batches, branch-free (also traced), the others select a component by the
# drawn threshold:
rand_impl(ctx::GenContext, μ::SuperpositionMeasure) = _superpose_rand(ctx, μ, _stream_ndims(μ))
_superpose_rand(ctx::GenContext, μ::SuperpositionMeasure, k::StaticInteger) = _superpose_batched_rand(ctx, μ, (), k)
function _superpose_rand(ctx::GenContext, μ::SuperpositionMeasure, ::NoMSpaceElementSize)
    components = values(μ.components)
    masses, total = _component_masses(μ)
    threshold = rand(ctx) * total
    csum = zero(threshold)
    for (mass, c) in zip(masses, components)
        csum += mass
        csum >= threshold && return rand_impl(ctx, c)
    end
    return rand_impl(ctx, last(components))
end

# Batches of superpositions draw a batch from each component and select
# by mass, branch-free. The draws are promoted to a common number type,
# point masses stay exact:
function batched_rand_impl(ctx::GenContext, μ::SuperpositionMeasure, sz::SizeLike)
    _superpose_batched_rand(ctx, μ, sz, _stream_ndims(μ))
end
function _superpose_batched_rand(ctx::GenContext, μ::SuperpositionMeasure, sz::SizeLike, k::StaticInteger)
    components = values(μ.components)
    masses, total = _component_masses(μ)
    thresholds = _batch_mask(_rand_bulk(ctx, sz) .* total, k)
    draws = map(c -> _component_draws(ctx, c, sz), components)
    T = mapreduce(eltype, promote_type, draws)
    Xs = map(X -> _convert_draws(T, X), draws)
    return _select_draws(thresholds, first(Xs), first(masses), _tail(masses), _tail(Xs))
end

@inline _component_draws(ctx::GenContext, μ, ::Tuple{}) = rand_impl(ctx, μ)
@inline _component_draws(ctx::GenContext, μ, sz::SizeLike) = batched_rand_impl(ctx, μ, sz)

@inline _convert_draws(::Type{T}, X::T) where {T<:Number} = X
@inline _convert_draws(::Type{T}, X::AbstractArray{T}) where {T} = X
@inline _convert_draws(::Type{T}, X) where {T} = convert.(T, X)

@inline _tail(xs::Tuple) = Base.tail(xs)
@inline _tail(xs::AbstractVector) = view(xs, firstindex(xs)+1:lastindex(xs))

# The first component whose cumulative mass reaches the threshold:
@inline _select_draws(thresholds, X, ::Any, ::Tuple{}, ::Tuple{}) = X
@inline function _select_draws(thresholds, X, csum, masses::Tuple, Ys::Tuple)
    X_new = ifelse.(thresholds .<= csum, X, first(Ys))
    _select_draws(thresholds, X_new, csum + first(masses), Base.tail(masses), Base.tail(Ys))
end
function _select_draws(thresholds, X, csum, masses::AbstractVector, Ys::AbstractVector)
    for (mass, Y) in zip(masses, Ys)
        X = ifelse.(thresholds .<= csum, X, Y)
        csum += mass
    end
    return X
end
_superpose_batched_rand(ctx::GenContext, μ::SuperpositionMeasure, sz::SizeLike, ::NoMSpaceElementSize) =
    _batched_rand_pointwise(ctx, μ, sz)

@inline function insupportof(d::SuperpositionMeasure, x)
    mapreduce(c -> _insupport_mask(insupportof(c, x)), |, values(d.components))
end


@inline mspace_flatsize(μ::SuperpositionMeasure) = mspace_flatsize(typeof(μ))

# The variate rank of a superposition is the common rank of its components,
# folded pairwise over the component types so that it stays a constant:
@inline mspace_ndims(::Type{<:SuperpositionMeasure{C}}) where {C<:AbstractArray} = mspace_ndims(eltype(C))
@inline function mspace_ndims(::Type{MU}) where {C<:Tuple,MU<:SuperpositionMeasure{C}}
    tuple_type_mapreduce(mspace_ndims, _CommonNDims{MU}(), C)
end
@inline stream_ndims(::Type{<:SuperpositionMeasure{C}}) where {C<:AbstractArray} = stream_ndims(eltype(C))
@inline function stream_ndims(::Type{MU}) where {C<:Tuple,MU<:SuperpositionMeasure{C}}
    tuple_type_mapreduce(stream_ndims, _CommonNDims{MU}(), C)
end
struct _CommonNDims{MU} <: Function end
@inline (::_CommonNDims{MU})(a::Integer, b::Integer) where {MU} = a == b ? a : NoMSpaceElementSize{MU}()
@inline (::_CommonNDims{MU})(::Any, ::Any) where {MU} = NoMSpaceElementSize{MU}()

@inline mspace_flatsize(::Type{MU}) where {C<:AbstractArray,MU<:SuperpositionMeasure{C}} = _scalar_or_unknown(MU, mspace_flatsize(eltype(C)))
@inline function mspace_flatsize(::Type{MU}) where {C<:Tuple,MU<:SuperpositionMeasure{C}}
    _scalar_or_unknown(MU, tuple_type_mapreduce(mspace_flatsize, _common_flatsize, C))
end
@inline _common_flatsize(a, b) = a === b ? a : NoMSpaceElementSize{typeof((a, b))}()
