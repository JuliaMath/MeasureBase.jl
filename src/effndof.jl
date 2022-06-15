"""
    effndof(μ)

Returns the effective number of degrees of freedom of variates of
measure `μ`.

The effective NDOF my differ from the length of the variates. For example,
the effective NDOF for a Dirichlet distribution with variates of length `n`
is `n - 1`.

Also see [`require_same_effndof`](@ref).
"""
function effndof end


"""
    MeasureBase.require_same_effndof(a, b)::Nothing

Check if `a` and `b` have the same effective number of degrees of freedom
according to [`MeasureBase.effndof`](@ref).
"""
function require_same_effndof end

ChainRulesCore.rrule(::typeof(require_same_effndof), a, b) = nothing, _nogradient_pullback2

function require_same_effndof(a, b)
    trg_d_n = effndof(ν)
    src_d_n = effndof(μ)
    if trg_d_n != src_d_n
        throw(ArgumentError("Can't convert to $(typeof(ν).name) with $(trg_d_n) eff. DOF from $(typeof(μ).name) with $(src_d_n) eff. DOF"))
    end
    return nothing
end
