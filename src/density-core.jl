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
"""
@inline logdensityof(μ::AbstractMeasure, x) = logdensityof_impl(μ, x)

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
    result = dynamic(unsafe_logdensityof(μ, x))
    _checksupport(insupport(μ, x), result)
end

@inline function logdensityof_rt(::T, ::U) where {T,U}
    Core.Compiler.return_type(logdensityof, Tuple{T,U})
end

_checksupport(cond, result) = ifelse(cond == true, result, oftype(result, -Inf))
@inline _checksupport(::NoFastInsupport, result) = result

"""
    MeasureBase.logdensityof_with_rest(μ::AbstractMeasure, x)

Compute the log-density of `μ` at the beginning of `x`, a flat stream of
variate content that may extend beyond the variate of `μ`.

`x` must either be a vector that starts with the (flattened) variate of `μ`
(for measures combined via `vcat`) or a `NamedTuple` whose first properties
constitute the variate of `μ` (for measures combined via `merge`).

Returns a tuple `(ℓ, x_μ, x_rest)` of the log-density `ℓ`, the variate
`x_μ` of `μ` consumed from the stream, and the unconsumed rest of the
stream.

Measure types whose variate size depends on measure values, like
[`mbind`](@ref) results, implement density calculation via
`logdensityof_with_rest` instead of
[`logdensityof_impl`](@ref MeasureBase.logdensityof_impl).

The default implementation determines the size resp. the property names of
the variate via [`some_mspace_elsize`](@ref MeasureBase.some_mspace_elsize)
resp. `testvalue` and delegates to `logdensityof_impl`.
"""
function logdensityof_with_rest end

function logdensityof_with_rest(μ::AbstractMeasure, x::AbstractVector)
    a, x_rest = _consume_from_stream(x, some_mspace_elsize(μ))
    return logdensityof_impl(μ, a), a, x_rest
end

function logdensityof_with_rest(μ::AbstractMeasure, x::NamedTuple)
    a, x_rest = _split_after(x, Val(_mspace_names(μ)))
    return logdensityof_impl(μ, a), a, x_rest
end

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
        ℓ_{i} = let Δℓ_{i} = logdensity_def(b_{i}, x)
            ℓ_{i - 1} + Δℓ_{i}
        end
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
@inline function logdensity_rel(μ::M, ν::N, x::X) where {M,N,X}
    inμ = insupport(μ, x)
    inν = insupport(ν, x)
    return _logdensity_rel_impl(μ, ν, x, inμ, inν)
end

@inline function _logdensity_rel_impl(μ::M, ν::N, x::X, inμ::Bool, inν::Bool) where {M,N,X}
    T = unstatic(
        promote_type(
            return_type(logdensity_def, (μ, x)),
            return_type(logdensity_def, (ν, x)),
        ),
    )
    istrue(inμ) || return convert(T, ifelse(inν, -Inf, NaN))
    istrue(inν) || return convert(T, Inf)

    return unsafe_logdensity_rel(μ, ν, x)
end

@inline function _logdensity_rel_impl(
    μ::M,
    ν::N,
    x::X,
    @nospecialize(::NoFastInsupport),
    @nospecialize(::NoFastInsupport)
) where {M,N,X}
    unsafe_logdensity_rel(μ, ν, x)
end

@inline function _logdensity_rel_impl(
    μ::M,
    ν::N,
    x::X,
    inμ::Bool,
    @nospecialize(::NoFastInsupport)
) where {M,N,X}
    logd = unsafe_logdensity_rel(μ, ν, x)
    return istrue(inμ) ? logd : oftype(logd, -Inf)
end

@inline function _logdensity_rel_impl(
    μ::M,
    ν::N,
    x::X,
    @nospecialize(::NoFastInsupport),
    inν::Bool
) where {M,N,X}
    logd = unsafe_logdensity_rel(μ, ν, x)
    return istrue(inν) ? logd : oftype(logd, +Inf)
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
