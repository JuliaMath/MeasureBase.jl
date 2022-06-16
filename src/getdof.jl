"""
    getdof(μ)

Returns the effective number of degrees of freedom of variates of
measure `μ`.

The effective NDOF my differ from the length of the variates. For example,
the effective NDOF for a Dirichlet distribution with variates of length `n`
is `n - 1`.

Also see [`check_dof`](@ref).
"""
function getdof end


"""
    MeasureBase.check_dof(a, b)::Nothing

Check if `a` and `b` have the same effective number of degrees of freedom
according to [`MeasureBase.getdof`](@ref).
"""
function check_dof end

ChainRulesCore.rrule(::typeof(check_dof), a, b) = nothing, _nogradient_pullback2

function check_dof(ν, μ)
    trg_d_n = getdof(ν)
    src_d_n = getdof(μ)
    if trg_d_n != src_d_n
        throw(ArgumentError("Can't convert to $(typeof(ν).name) with $(trg_d_n) eff. DOF from $(typeof(μ).name) with $(src_d_n) eff. DOF"))
    end
    return nothing
end
