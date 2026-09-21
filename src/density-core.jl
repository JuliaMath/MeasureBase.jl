export logdensityof
export logdensity_rel
export logdensity_def

export unsafe_logdensityof
export unsafe_logdensity_rel

export densityof
export density_rel
export density_def

"""
    logdensityof(m::AbstractMeasure, x)

Compute the log-density of the measure `m` at `x`. Density is always relative,
but `DensityInterface.jl` does not account for this. For compatibility with
this, `logdensityof` for a measure is always implicitly relative to
[`rootmeasure(x)`](@ref rootmeasure).

`logdensityof(m, x)` is implemented via
[`MeasureBase.logdensityof_impl`](@ref), measure types should specialize
`logdensityof_impl` instead of `logdensityof` itself.

To compute log-density relative to `basemeasure(m)` or *define* a log-density
(relative to `basemeasure(m)` or another measure given explicitly), see
`logdensity_def`.

To compute a log-density relative to a specific base-measure, see
`logdensity_rel`.

# Extended help

Variates of the right shape and element type never throw: outside the
support of `m` the result is `-Inf`, also for non-integer values of
measures over counting measures and for infinite values. `NaN` inputs give
`NaN` or `-Inf`. Variates of the wrong shape throw an `ArgumentError`.
Implementations of `logdensityof_impl` and `unsafe_logdensityof` must not
throw outside the support, since support masks evaluate both branches.
"""
@inline logdensityof(μ::AbstractMeasure, x) = _point_ld(logdensityof_impl, μ, x)

"""
    MeasureBase.logdensityof_impl(μ::AbstractMeasure, x)

Implements [`logdensityof(μ, x)`](@ref logdensityof).

Measure types should specialize `logdensityof_impl` instead of
`logdensityof` itself. Implementations must return the log-density of `μ`
at `x` relative to [`rootmeasure(μ)`](@ref) and must handle `x` outside of
the support of `μ` (the result must be `-Inf` then).

The default implementation checks `insupport(μ, x)` (unless the result is
a [`MeasureBase.NoFastInsupport`](@ref)) and computes the density via
[`unsafe_logdensityof`](@ref).
"""
@inline function logdensityof_impl(μ::AbstractMeasure, x)
    result = _dynamic_logd(unsafe_logdensityof(μ, x), x)
    _checksupport(insupport(μ, x), result)
end

# Log-density kernels return numbers of the number type of the variate,
# never static numbers, so that automatic differentiation and tracing see
# ordinary floating point values throughout:
@inline _logd_numtype(x) = float(real_numtype(typeof(x)))
@inline _dynamic_logd(ℓ, x) = dynamic(ℓ)
@inline _neg_inf_logd(x) = _logd_numtype(x)(-Inf)

# Support checks as masks: `NoFastInsupport` means the density is evaluated
# unconditionally.
@inline _insupport_mask(ins) = ins == true
@inline _insupport_mask(::NoFastInsupport) = true

# Support checks as booleans, keeping `NoFastInsupport`:
@inline _insupport_bool(ins) = ins == true
@inline _insupport_bool(ins::NoFastInsupport) = ins

# Combining support checks of components, `NoFastInsupport` is absorbing:
@inline _insupport_and(a, b) = _insupport_bool(a) & _insupport_bool(b)
@inline _insupport_and(a::NoFastInsupport, ::Any) = a
@inline _insupport_and(::Any, b::NoFastInsupport) = b
@inline _insupport_and(a::NoFastInsupport, ::NoFastInsupport) = a

@inline _checksupport(cond, result) = ifelse(_insupport_mask(cond), result, oftype(result, -Inf))

# Transports of variates outside the support of the source measure give
# NaN. Both branches are evaluated, formulas must not throw outside the
# support:
@inline _nan_outside(μ, x, y) = ifelse(_insupport_mask(insupport(μ, x)), y, oftype(y, NaN))

# On the floating-point grid the endpoints of the unit interval stand for
# their nearest interior points (the smallest normal float above zero,
# since devices may flush subnormals, and the grid point below one), so
# that quantiles stay finite, and tail probabilities in log-space
# conversions never underflow to zero:
@inline _unit_interior(p) = clamp(p, _unit_bounds(p)...)
@inline _unit_bounds(p) = (_prob_floor(p), prevfloat(one(p)))
@inline _positive_prob(p) = max(p, _prob_floor(p))
@inline _prob_floor(p) = floatmin(typeof(one(p)))

"""
    MeasureBase.logdensityof_with_rest(μ::AbstractMeasure, x)

Consume the variate of `μ` at the beginning of the flat vector stream `x`
(or the named entries of the `NamedTuple` `x`) and compute its log-density.

Returns a tuple `(ℓ, x_μ, x_rest)` of the log-density, the consumed variate
`x_μ` and the unconsumed rest of `x`. Measures whose variates have a fixed
size consume that size (see [`MeasureBase.mspace_flatsize`](@ref) and
[`MeasureBase.some_mspace_elsize`](@ref)), measures with variates of
value-dependent size implement the consumption themselves. Batches of
streams are consumed by [`MeasureBase.batched_logdensityof_with_rest`](@ref).
"""
function logdensityof_with_rest end

function logdensityof_with_rest(μ::AbstractMeasure, x::AbstractVector)
    a, x_rest = _consume_from_stream(x, _stream_consume_size(μ))
    return _point_ld(logdensityof_impl, μ, a), a, x_rest
end

function logdensityof_with_rest(μ::AbstractMeasure, x::NamedTuple)
    a, x_rest = _split_after(x, Val(_mspace_names(μ)))
    return logdensityof_impl(μ, a), a, x_rest
end

@inline _stream_consume_size(μ) = _stream_consume_size(μ, mspace_flatsize(μ))
@inline _stream_consume_size(μ, sz::SizeLike) = sz
@inline _stream_consume_size(μ, ::NoMSpaceElementSize) = some_mspace_elsize(μ)

_mspace_names(μ::AbstractMeasure) = keys(testvalue(μ))


"""
    localmeasure(m::AbstractMeasure, x)::AbstractMeasure

Return a measure that behaves like `m` in the infinitesimal neighborhood
of `x` in respect to density calculation.

Note that the resulting measure may not be well defined outside of the
infinitesimal neighborhood of `x`.

For most measure types simply returns `m` itself. [`mbind`](@ref),
for example, generates measures for which `localmeasure(m, x)` depends
on `x`.
"""
localmeasure(m::AbstractMeasure, x) = m
export localmeasure


"""
    MeasureBase.transportmeasure(m::AbstractMeasure, x)::AbstractMeasure

Return a measure that behaves like `m` in the infinitesimal neighborhood
of `x` in respect to both transport and density calculation.

Note that the resulting measure may not be well defined outside of the
infinitesimal neighborhood of `x`.

For most measure types simply returns `m` itself. [`mbind`](@ref),
for example, generates measures for which `transportmeasure(m, x)` depends
on `x`.
"""
transportmeasure(m::AbstractMeasure, x) = m

export unsafe_logdensityof

# https://discourse.julialang.org/t/counting-iterations-to-a-type-fixpoint/75876/10?u=cscherrer
"""
    unsafe_logdensityof(m, x)

Compute the log-density of the measure `m` at `x` relative to `rootmeasure(m)`.
This is "unsafe" because it does not check `insupport(m, x)`.

See also `logdensityof`.
"""
@inline function unsafe_logdensityof(μ::AbstractMeasure, x)
    μ_local = localmeasure(μ, x)
    # Extra dispatch boundary to reduce number of required specializations of implementation:
    return _unsafe_logdensityof_local(μ_local, x)
end

@inline function _unsafe_logdensityof_local(μ::M, x) where {M}
    ℓ_0 = logdensity_def(μ, x)
    b_0 = μ
    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        b_{i} = basemeasure(b_{i - 1})

        # The below makes the evaluated code shorter, but screws up Zygote
        # if b_{i} isa typeof(b_{i - 1})
        #     return ℓ_{i - 1}
        # end
        ℓ_{i} = ℓ_{i - 1} + logdensity_def(b_{i}, x)
    end
    return ℓ_10
end

"""
    logdensity_rel(m1, m2, x)

Compute the log-density of `m1` relative to `m2` at `x`. This function checks
whether `x` is in the support of `m1` or `m2` (or both, or neither). If `x` is
known to be in the support of both, it can be more efficient to call
`unsafe_logdensity_rel`. 
"""
@inline function logdensity_rel(μ, ν, x)
    inμ = _insupport_mask(insupport(μ, x))
    inν = _insupport_mask(insupport(ν, x))
    logd = _dynamic_logd(unsafe_logdensity_rel(μ, ν, x), x)
    outside = ifelse(inμ, oftype(logd, +Inf), ifelse(inν, oftype(logd, -Inf), oftype(logd, NaN)))
    return ifelse(inμ & inν, logd, outside)
end

"""
    unsafe_logdensity_rel(m1, m2, x)

Compute the log-density of `m1` relative to `m2` at `x`, assuming `x` is
known to be in the support of both `m1` and `m2`.

See also `logdensity_rel`.
"""
@inline function unsafe_logdensity_rel(μ::AbstractMeasure, ν::AbstractMeasure, x)
    μ_local = localmeasure(μ, x)
    ν_local = localmeasure(ν, x)
    return logdensity_def(μ_local, ν_local, x)
end

# Indicates that no specialized method is available to compute the
# log-density between a given pair of measures:
struct _NoLogdensityRel end

"""
    MeasureBase.logdensity_rel_def(μ, ν, x)

Specialization point for the log-density of `μ` relative to `ν` at `x`.

Measure types may add methods for pairs of measure types whose relative
density can be computed directly. The generic implementation of
[`logdensity_def(μ, ν, x)`](@ref logdensity_def) descends the base measure
chains of both measures in lockstep and uses the first specialized
`logdensity_rel_def` method it encounters along the way.

Do not call `logdensity_rel_def` directly, call
[`logdensity_rel`](@ref) (or `logdensity_def`) instead.
"""
@inline logdensity_rel_def(μ, ν, x) = _NoLogdensityRel()

# Generic relative density: descend the base measure chains of both measures
# in lockstep, after equalizing their depths. Since the members of a shared
# chain suffix have the same depth-from-root on both sides, the descent
# terminates at a specialized `logdensity_rel_def` method as soon as one
# becomes applicable (in particular for pairs of identical primitive
# measures), so any shared chain suffix cancels symbolically instead of
# numerically. The descent is fully unrolled at compile time based on the
# static base measure depths, only the base measures actually visited are
# constructed, and whether a specialized method applies at a given level is
# decided purely by dispatch (on the sentinel type `_NoLogdensityRel`).
@inline function logdensity_def(μ, ν, x)
    _logdensity_rel_descent(μ, basemeasure_depth(μ), ν, basemeasure_depth(ν), x)
end

@generated function _logdensity_rel_descent(
    μ,
    ::StaticInteger{M},
    ν,
    ::StaticInteger{N},
    x,
) where {M,N}
    μsym(i) = Symbol(:μ_, i)
    νsym(j) = Symbol(:ν_, j)
    prog = Expr(:block, Expr(:meta, :inline), :(μ_0 = μ), :(ν_0 = ν))
    terms = Any[]
    n_checks = 0
    # Return via a specialized `logdensity_rel_def` method for the current
    # measure pair, if available. Whether one is available is decided purely
    # by type, so unsuccessful checks are free at run time:
    function emit_check!(i, j)
        r = Symbol(:r_, n_checks)
        n_checks += 1
        push!(prog.args, :($r = logdensity_rel_def($(μsym(i)), $(νsym(j)), x)))
        ret = isempty(terms) ? r : :(+($(terms...), $r))
        push!(prog.args, :(if !($r isa _NoLogdensityRel)
            return $ret
        end))
    end
    i = j = 0
    emit_check!(i, j)
    # Equalize depths, accumulating one-sided density terms:
    while M - i > N - j
        ℓ = Symbol(:ℓμ_, i)
        push!(prog.args, :($ℓ = logdensity_def($(μsym(i)), x)))
        push!(prog.args, :($(μsym(i + 1)) = basemeasure($(μsym(i)))))
        push!(terms, ℓ)
        i += 1
        emit_check!(i, j)
    end
    while N - j > M - i
        ℓ = Symbol(:ℓν_, j)
        push!(prog.args, :($ℓ = -logdensity_def($(νsym(j)), x)))
        push!(prog.args, :($(νsym(j + 1)) = basemeasure($(νsym(j)))))
        push!(terms, ℓ)
        j += 1
        emit_check!(i, j)
    end
    # Lockstep descent at equal depth:
    for _ in 1:(M-i)
        ℓμ, ℓν = Symbol(:ℓμ_, i), Symbol(:ℓν_, j)
        push!(prog.args, :($ℓμ = logdensity_def($(μsym(i)), x)))
        push!(prog.args, :($ℓν = -logdensity_def($(νsym(j)), x)))
        push!(terms, ℓμ, ℓν)
        push!(prog.args, :($(μsym(i + 1)) = basemeasure($(μsym(i)))))
        push!(prog.args, :($(νsym(j + 1)) = basemeasure($(νsym(j)))))
        i += 1
        j += 1
        emit_check!(i, j)
    end
    # Both measures are at root level now:
    push!(
        prog.args,
        :(r_root = _root_logdensity_rel($(μsym(i)), $(νsym(j)), x)),
    )
    ret = isempty(terms) ? :r_root : :(+($(terms...), r_root))
    push!(prog.args, :(return $ret))
    return prog
end

# Root measures of the same type are equal almost everywhere for the
# purpose of pointwise relative densities:
_root_logdensity_rel(μ::M, ν::M, x) where {M} = zero(logdensity_def(μ, x))

function _root_logdensity_rel(@nospecialize(μ), @nospecialize(ν), @nospecialize(x))
    throw(
        ArgumentError(
            "No method available to compute the log-density between measures with root measures of type $(nameof(typeof(μ))) and $(nameof(typeof(ν)))",
        ),
    )
end

@inline density_rel(μ, ν, x) = exp(logdensity_rel(μ, ν, x))

# TODO: Do we need this method?
density_def(μ, ν::AbstractMeasure, x) = exp(logdensity_def(μ, ν, x))
density_def(μ, x) = exp(logdensity_def(μ, x))
