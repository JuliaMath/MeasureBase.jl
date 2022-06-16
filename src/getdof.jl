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
    MeasureBase.check_dof(ν, μ)::Nothing

Check if `ν` and `μ` have the same effective number of degrees of freedom
according to [`MeasureBase.getdof`](@ref).
"""
function check_dof end

function check_dof(ν, μ)
    n_ν = getdof(ν)
    n_μ = getdof(μ)
    if n_ν != n_μ
        throw(ArgumentError("Measure ν of type $(nameof(typeof(ν))) has $(n_ν) DOF but μ of type $(nameof(typeof(μ))) has $(n_μ) DOF"))
    end
    return nothing
end

ChainRulesCore.rrule(::typeof(check_dof), ν, μ) = NoTangent(), NoTangent(), NoTangent()


"""
    MeasureBase.NoVarCheck{MU,T}

Indicates that there is no way to check of a values of type `T` are
variate of measures of type `MU`.
"""
struct NoVarCheck{MU,T} end


"""
    MeasureBase.checked_var(μ::MU, x::T)::T

Return `x` if `x` is a valid variate of `μ`, throw an `ArgumentError` if not,
return `NoVarCheck{MU,T}()` if not check can be performed.
"""
function checked_var end

@inline checked_var(::MU, ::T) where {MU,T} = NoVarCheck{MU,T}

ChainRulesCore.rrule(::typeof(checked_var), ν, x) = NoTangent(), NoTangent(), ZeroTangent()
