"""
    MeasureBase.NoFastInsupport{MU}

Indicates that there is no fast way to compute if a point lies within the
support of measures of type `MU`.
"""
struct NoFastInsupport{MU} end


"""
    insupportof(m, x)
    insupportof(m)

`insupportof(m, x)` computes whether `x` is in the support of `m` and
returns either a `Bool` (possibly a static one) or an instance of
[`NoFastInsupport`](@ref). Densities are `-Inf` outside of the support.

`insupportof(m)` returns a function, and satisfies
`insupportof(m)(x) == insupportof(m, x)`.
"""
function insupportof end

@deprecate insupport(m, x) insupportof(m, x)
@deprecate insupport(m) insupportof(m)

"""
    MeasureBase.require_insupport(μ, x)::Nothing

Checks if `x` is in the support of distribution/measure `μ`, throws an
`ArgumentError` if not.

Will not throw an exception if `insupportof` returns an instance of
[`NoFastInsupport`](@ref).
"""
function require_insupport end

function require_insupport(μ, x)
    ins = insupportof(μ, x)
    if !(ins isa NoFastInsupport)
        if !ins
            throw(ArgumentError("x is not within the support of μ"))
        end
    end
    return nothing
end
