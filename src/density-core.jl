export logdensityof
export logdensity_rel
export logdensity_def

export unsafe_logdensityof
export unsafe_logdensity_rel

export densityof
export density_rel
export density_def


"""
    localmeasure(m::AbstractMeasure, x)::AbstractMeasure

Return a measure that behaves like `m` in the infinitesimal neighborhood
of `x` in respect to density calculation.

Note that the resulting measure may not be well defined outside of the
infinitesimal neighborhood of `x`.

For most measure types simply returns `m` itself. [`mbind`](@ref),
for example, generates measures for with `localmeasure(m, x)` depends
on `x`.
"""
localmeasure(m::AbstractMeasure, x) = m
export localmeasure


"""
    MeasureBase.transportmeasure(μ::Bind, x)::AbstractMeasure

Return a measure that behaves like `m` in the infinitesimal neighborhood
of `x` in respect to both transport and density calculation.

Note that the resulting measure may not be well defined outside of the
infinitesimal neighborhood of `x`.

For most measure types simply returns `m` itself. [`mbind`](@ref),
for example, generates measures for with `transportmeasure(m, x)` depends
on `x`.
"""
transportmeasure(m::AbstractMeasure, x) = m
export localmeasure


"""
    logdensityof(m::AbstractMeasure, x) 

Compute the log-density of the measure `m` at `x`. Density is always relative,
but `DensityInterface.jl` does not account for this. For compatibility with
this, `logdensityof` for a measure is always implicitly relative to
[`rootmeasure(x)`](@ref rootmeasure). 

`logdensityof` works by first computing `insupport(m, x)`. If this is true, then
`unsafe_logdensityof` is called. If `insupport(m, x)` is known to be `true`, it
can be a little faster to directly call `unsafe_logdensityof(m, x)`. 

To compute log-density relative to `basemeasure(m)` or *define* a log-density
(relative to `basemeasure(m)` or another measure given explicitly), see
`logdensity_def`. 

To compute a log-density relative to a specific base-measure, see
`logdensity_rel`. 
"""
@inline function logdensityof(μ::AbstractMeasure, x)
    result = dynamic(unsafe_logdensityof(μ, x))
    _checksupport(insupport(μ, x), result)
end

_checksupport(cond, result) = ifelse(cond == true, result, oftype(result, -Inf))
@inline _checksupport(::NoFastInsupport, result) = result


export unsafe_logdensityof

# https://discourse.julialang.org/t/counting-iterations-to-a-type-fixpoint/75876/10?u=cscherrer
"""
    unsafe_logdensityof(m, x)

Compute the log-density of the measure `m` at `x` relative to `rootmeasure(m)`.
This is "unsafe" because it does not check `insupport(m, x)`.

See also `logdensityof`.
"""
@inline function unsafe_logdensityof(μ::M, x) where {M}
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
    logdensity_type(m::AbstractMeasure}, ::Type{T}) where T

Compute the return type of `logdensity_of(m, ::T)`.
"""
function logdensity_type(m::M,T) where {M<:AbstractMeasure}
    Core.Compiler.return_type(logdensity_def, Tuple{M, T})
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
            logdensity_type(μ, X),
            logdensity_type(ν, X),
        ),
    )

    istrue(inμ) || return convert(T, ifelse(inν, -Inf, NaN))
    istrue(inν) || return convert(T, Inf)

    return unsafe_logdensity_rel(μ, ν, x)
end


@inline function _logdensity_rel_impl(μ::M, ν::N, x::X, @nospecialize(::NoFastInsupport), @nospecialize(::NoFastInsupport)) where {M,N,X}
    unsafe_logdensity_rel(μ, ν, x)
end

@inline function _logdensity_rel_impl(μ::M, ν::N, x::X, inμ::Bool, @nospecialize(::NoFastInsupport)) where {M,N,X}
    logd = unsafe_logdensity_rel(μ, ν, x)
    return istrue(inμ) ? logd  : logd * oftype(logd, -Inf)
end

@inline function _logdensity_rel_impl(μ::M, ν::N, x::X, @nospecialize(::NoFastInsupport), inν::Bool) where {M,N,X}
    logd = unsafe_logdensity_rel(μ, ν, x)
    return istrue(inν) ? logd  : logd * oftype(logd, +Inf)
end


"""
    unsafe_logdensity_rel(m1, m2, x)

Compute the log-density of `m1` relative to `m2` at `x`, assuming `x` is
known to be in the support of both `m1` and `m2`.

See also `logdensity_rel`.
"""
@inline function unsafe_logdensity_rel(μ::M, ν::N, x::X) where {M,N,X}
    μ_local = localmeasure(μ, x)
    ν_local = localmeasure(ν, x)
    # Extra dispatch boundary to reduce number of required specializations of implementation:
    return _unsafe_logdensity_rel_local(μ_local, ν_local, x)
end

@inline function _unsafe_logdensity_rel_local(μ::M, ν::N, x::X) where {M,N,X}
    if static_hasmethod(logdensity_def, Tuple{M,N,X})
        return logdensity_def(μ, ν, x)
    end
    μs = basemeasure_sequence(μ)
    νs = basemeasure_sequence(ν)
    cb = commonbase(μs, νs, X)
    # _logdensity_rel(μ, ν)
    isnothing(cb) && begin
        μ = μs[end]
        ν = νs[end]
        @warn """
        No common base measure for
            $μ
        and
            $ν

        Returning a relative log-density of NaN. If this is incorrect, add a
        three-argument method
            logdensity_def($μ, $ν, x)
        """
        return NaN
    end
    return _logdensity_rel(μs, νs, cb, x)
end

# Note that this method assumes `μ` and `ν` to have the same type
function logdensity_def(μ::T, ν::T, x) where {T}
    if μ === ν
        return zero(logdensity_def(μ, x))
    else
        α = basemeasure(μ)
        β = basemeasure(ν)
        return logdensity_def(μ, x) - logdensity_def(ν, x) + logdensity_rel(α, β, x)
    end
end

@generated function _logdensity_rel(
    μs::Tμ,
    νs::Tν,
    ::Tuple{<:StaticInteger{M},<:StaticInteger{N}},
    x::X,
) where {Tμ,Tν,M,N,X}
    sμ = schema(Tμ)
    sν = schema(Tν)

    q = quote
        $(Expr(:meta, :inline))
        ℓ = logdensity_def(μs[$M], νs[$N], x)
    end

    for i in 1:M-1
        push!(q.args, :(Δℓ = logdensity_def(μs[$i], x)))
        # push!(q.args, :(println("Adding", Δℓ)))
        push!(q.args, :(ℓ += Δℓ))
    end

    for j in 1:N-1
        push!(q.args, :(Δℓ = logdensity_def(νs[$j], x)))
        # push!(q.args, :(println("Subtracting", Δℓ)))
        push!(q.args, :(ℓ -= Δℓ))
    end

    push!(q.args, :(return ℓ))
    return q
end

@inline density_rel(μ, ν, x) = exp(logdensity_rel(μ, ν, x))

# TODO: Do we need this method?
density_def(μ, ν::AbstractMeasure, x) = exp(logdensity_def(μ, ν, x))
density_def(μ, x) = exp(logdensity_def(μ, x))
