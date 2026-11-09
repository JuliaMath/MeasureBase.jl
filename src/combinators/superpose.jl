
using LogarithmicNumbers
using LogExpFunctions

export SuperpositionMeasure

abstract type AbstractSuperpositionMeasure <: AbstractMeasure end

"""
    MeasureBase.MeasureJoin(roots)

The join of the root measures `roots` (a tuple or an array of measures),
the smallest measure that dominates all of them. Its density relative to
each of the roots is the largest density of the roots relative to that
one: equal roots join to the root itself, mutually singular roots (like
Lebesgue measure and a point mass) to their sum.

`MeasureJoin` is the base and root measure of superpositions (see
[`SuperpositionMeasure`](@ref)) and relates to other measures through the
roots it joins. Joins of arrays of roots record at construction whether
all roots are equal.
"""
struct MeasureJoin{R} <: AbstractMeasure
    roots::R
    equal_roots::Bool

    MeasureJoin(roots::Tuple) = new{typeof(roots)}(roots, false)
    MeasureJoin(roots::AbstractArray) = new{typeof(roots)}(roots, allequal(roots))
end

Base.:(==)(a::MeasureJoin, b::MeasureJoin) = a.roots == b.roots
Base.isequal(a::MeasureJoin, b::MeasureJoin) = isequal(a.roots, b.roots)
Base.hash(a::MeasureJoin, h::UInt) = hash(a.roots, hash(:MeasureJoin, h))

function Pretty.tile(J::MeasureJoin)
    Pretty.literal("MeasureJoin(") * Pretty.list_layout([Pretty.tile.(J.roots)...]) * Pretty.literal(")")
end

# Roots of identical singleton types are equal, their join is the root:
@inline _measure_join(roots::Tuple{R,Vararg{R}}) where {R} = Base.issingletontype(R) ? first(roots) : MeasureJoin(roots)
@inline _measure_join(roots::Tuple) = MeasureJoin(roots)
@inline _measure_join(roots::AbstractArray{R}) where {R} = Base.issingletontype(R) ? instance(R) : MeasureJoin(roots)

basemeasure(J::MeasureJoin) = J
logdensity_def(::MeasureJoin, x) = zero(_logd_numtype(x))

@inline insupportof(J::MeasureJoin, x) = mapreduce(r -> _insupport_mask(insupportof(r, x)), |, J.roots)

# The log-density of the join relative to a measure ν, the largest
# log-density of the joined roots relative to ν, roots without mass at `x`
# excluded:
@inline function _log_djoin(J::MeasureJoin, ν, x)
    J.equal_roots && return logdensity_rel(J.roots[begin], ν, x)
    _max_logd(r -> _masked_logd(logdensity_rel(r, ν, x), insupportof(r, x)), J.roots)
end

@inline _max_logd(f::F, roots::Tuple) where {F} = reduce(max, map(f, roots))
@inline _max_logd(f::F, roots::AbstractArray) where {F} = mapreduce(f, max, roots)

@inline logdensity_rel_def(J::MeasureJoin, ν, x) = _log_djoin(J, ν, x)
@inline logdensity_rel_def(ν, J::MeasureJoin, x) = -_log_djoin(J, ν, x)
@inline function logdensity_rel_def(J::MeasureJoin, K::MeasureJoin, x)
    _max_logd(r -> _masked_logd(-_log_djoin(K, r, x), insupportof(r, x)), J.roots)
end

# Joins of roots of one kind have that kind (see `_root_kind`):
@inline _root_kind(J::MeasureJoin) = _common_root_kind(J.roots)
@inline _lebesgue_logd(J::MeasureJoin, x) = _max_logd(r -> _lebesgue_logd(r, x), J.roots)
@inline _atom_free(J::MeasureJoin) = all(_atom_free, J.roots)
@inline _pushfwd_root(J::MeasureJoin) = _any_pushfwd_root(J.roots)

@doc raw"""
    struct SuperpositionMeasure{NT} <: AbstractMeasure
        components :: NT
    end

Superposition of measures is analogous to mixture distributions, but (because
measures need not be normalized) requires no scaling.
The superposition of two measures μ and ν can be more concisely written as μ + ν.

The base and root measure of a superposition is the join of the root
measures ``r_i`` of its components (see [`MeasureBase.MeasureJoin`](@ref)).
Relative to it, the superposition has the density

```math
    \frac{\mathrm{d}\sum_i \mu_i}{\mathrm{d}\bigvee_j r_j}(x) = \sum_i
    \frac{\mathrm{d}\mu_i}{\mathrm{d}r_i}(x) \Big/ \max_j \frac{\mathrm{d}r_j}{\mathrm{d}r_i}(x),
```

the sum of the component densities for components with equal roots, the
weight of an atom at the atom for a point mass added to a continuous
measure.
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

# The density of a superposition relative to the join of the component
# roots, in log space: each component contributes its own density, divided
# by the density of the join relative to its own root. Components with
# equal roots add their densities:
logdensity_def(s::SuperpositionMeasure, x) = _superpos_logd(values(s.components), basemeasure(s), x)

_superpos_logd(cs, ::Any, x) = _logsumexp_components(map(c -> _component_logd(c, x), cs))
function _superpos_logd(cs, J::MeasureJoin, x)
    J.equal_roots && return _superpos_logd(cs, J.roots[begin], x)
    _logsumexp_components(map(cs, J.roots) do cᵢ, rᵢ
        _rel_to_join(_component_logd(cᵢ, x), _log_djoin(J, rᵢ, x))
    end)
end

@inline _component_logd(c, x) = _dynamic_logd(logdensityof(c, x), x)

# Components without mass at `x` contribute nothing, also where the join
# has no density relative to their root:
@inline _rel_to_join(ℓ, log_dJ) = ifelse(isneginf(ℓ), ℓ, ℓ - log_dJ)

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

@inline logdensity_rel_def(s::SuperpositionMeasure, J::MeasureJoin, x) = _superpos_logdensity_rel(s, J, x)
@inline logdensity_rel_def(J::MeasureJoin, s::SuperpositionMeasure, x) = -_superpos_logdensity_rel(s, J, x)

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

basemeasure(μ::SuperpositionMeasure) = _measure_join(map(rootmeasure, values(μ.components)))

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
@inline _convert_draws(::Type{T}, x::Number) where {T} = convert(T, x)
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
@inline mspace_flatsize(J::MeasureJoin) = mspace_flatsize(typeof(J))

testvalue(::Type{T}, J::MeasureJoin) where {T} = testvalue(T, J.roots[begin])

# The variate rank of a superposition (or a join) is the common rank of its
# components, folded pairwise over the component types so that it stays a
# constant:
const _MeasureCollection{C} = Union{SuperpositionMeasure{C},MeasureJoin{C}}

@inline mspace_ndims(::Type{<:_MeasureCollection{C}}) where {C<:AbstractArray} = mspace_ndims(eltype(C))
@inline function mspace_ndims(::Type{MU}) where {C<:Tuple,MU<:_MeasureCollection{C}}
    tuple_type_mapreduce(mspace_ndims, _CommonNDims{MU}(), C)
end
@inline stream_ndims(::Type{<:_MeasureCollection{C}}) where {C<:AbstractArray} = stream_ndims(eltype(C))
@inline function stream_ndims(::Type{MU}) where {C<:Tuple,MU<:_MeasureCollection{C}}
    tuple_type_mapreduce(stream_ndims, _CommonNDims{MU}(), C)
end
struct _CommonNDims{MU} <: Function end
@inline (::_CommonNDims{MU})(a::Integer, b::Integer) where {MU} = a == b ? a : NoMSpaceElementSize{MU}()
@inline (::_CommonNDims{MU})(::Any, ::Any) where {MU} = NoMSpaceElementSize{MU}()

@inline mspace_flatsize(::Type{MU}) where {C<:AbstractArray,MU<:_MeasureCollection{C}} = _scalar_or_unknown(MU, mspace_flatsize(eltype(C)))
@inline function mspace_flatsize(::Type{MU}) where {C<:Tuple,MU<:_MeasureCollection{C}}
    _scalar_or_unknown(MU, tuple_type_mapreduce(mspace_flatsize, _common_flatsize, C))
end
@inline _common_flatsize(a, b) = a === b ? a : NoMSpaceElementSize{typeof((a, b))}()
