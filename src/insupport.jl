"""
    inssupport(m, x)
    insupport(m)

`insupport(m,x)` computes whether `x` is in the support of `m`.

`insupport(m)` returns a function, and satisfies

insupport(m)(x) == insupport(m, x)
"""
function insupport end


"""
    MeasureBase.require_insupport(μ, x)::Nothing

Checks if `x` is in the support of distribution/measure `μ`, throws an
`ArgumentError` if not.
"""
function require_insupport end

_check_insupport_pullback(ΔΩ) = NoTangent(), ZeroTangent()
function ChainRulesCore.rrule(::typeof(require_insupport), μ, x)
    return require_insupport(μ, x), _check_insupport_pullback
end

function require_insupport(μ, x::AbstractArray{T,N}) where {T,N}
    if !insupport(μ, x)
        throw(ArgumentError("x is not within the support of μ"))
    end
    return nothing
end


"""
    MeasureBase.check_varshape(μ, x)::Nothing

Checks if `x` has the correct shape/size for a variate of measure-like object
`μ`, throws an `ArgumentError` if not.
"""
function check_varshape end

_check_varshape_pullback(ΔΩ) = NoTangent(), ZeroTangent()
ChainRulesCore.rrule(::typeof(check_varshape), μ, x) = check_varshape(μ, x), _check_varshape_pullback
