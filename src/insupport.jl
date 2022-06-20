"""
    inssupport(m, x)
    insupport(m)

`insupport(m,x)` computes whether `x` is in the support of `m`.

`insupport(m)` returns a function, and satisfies

insupport(m)(x) == insupport(m, x)
"""
function insupport end

"""
    require_insupport(μ, x)::Nothing

Checks if `x` is in the support of distribution/measure `μ`, throws an
`ArgumentError` if not.
"""
function require_insupport end

_require_insupport_pullback(ΔΩ) = NoTangent(), ZeroTangent()
function ChainRulesCore.rrule(::typeof(require_insupport), μ, x)
    return require_insupport(μ, x), _require_insupport_pullback
end

function require_insupport(μ, x)
    if !insupport(μ, x)
        throw(ArgumentError("x is not within the support of μ"))
    end
    return nothing
end
