"""
    MeasureBase.NoDOF{MU}

Indicates that there is no way to compute degrees of freedom of a measure
of type `MU` with the given information, e.g. because the DOF are not
a global property of the measure.
"""
struct NoDOF{MU} end

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

# Prevent infinite recursion:
@inline _default_getdof(::Type{MU}, ::MU) where {MU} = NoDOF{MU}
@inline _default_getdof(::Type{MU}, mu_base) where {MU} = getdof(mu_base)

@inline getdof(μ::MU) where {MU} = _default_getdof(MU, basemeasure(μ))

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
        throw(
            ArgumentError(
                "Measure ν of type $(nameof(typeof(ν))) has $(n_ν) DOF but μ of type $(nameof(typeof(μ))) has $(n_μ) DOF",
            ),
        )
    end
    return nothing
end

"""
    MeasureBase.NoArgCheck{MU,T}

Indicates that there is no way to check of a values of type `T` are
variate of measures of type `MU`.
"""
struct NoArgCheck{MU,T} end

"""
    MeasureBase.checked_arg(μ::MU, x::T)::T

Return `x` if `x` is a valid variate of `μ`, throw an `ArgumentError` if not,
return `NoArgCheck{MU,T}()` if not check can be performed.
"""
function checked_arg end

# Prevent infinite recursion:
@propagate_inbounds function _default_checked_arg(::Type{MU}, ::MU, ::T) where {MU,T}
    NoArgCheck{MU,T}
end
@propagate_inbounds function _default_checked_arg(::Type{MU}, mu_base, x) where {MU}
    checked_arg(mu_base, x)
end

@propagate_inbounds function checked_arg(mu::MU, x) where {MU}
    _default_checked_arg(MU, basemeasure(mu), x)
end
