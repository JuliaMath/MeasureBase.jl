"""
    MeasureBase.NoFastInsupport{MU}

Indicates that there is no fast way to compute if a point lies within the
support of measures of type `MU`
"""
struct NoFastInsupport{MU} end


"""
    inssupport(m, x)
    insupport(m)

`insupport(m,x)` computes whether `x` is in the support of `m` and
returns either a `Bool` or an instance of [`NoFastInsupport`](@ref).

`insupport(m)` returns a function, and satisfies
`insupport(m)(x) == insupport(m, x)``
"""
function insupport end


"""
    MeasureBase.require_insupport(μ, x)::Nothing

Checks if `x` is in the support of distribution/measure `μ`, throws an
`ArgumentError` if not.

Will not throw an exception if `insupport` returns an instance of
[`NoFastInsupport`](@ref).
"""
function require_insupport end

function require_insupport(μ, x)
    ins = insupport(μ, x)
    if !(ins isa NoFastInsupport)
        if !ins
            throw(ArgumentError("x is not within the support of μ"))
        end
    end
    return nothing
end
