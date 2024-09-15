"""
    abstract type MeasureBase.AbstractNoDOF{MU}

Abstract supertype for [`NoDOF`](@ref) and [`NoFastDOF`](@ref).
"""
abstract type AbstractNoDOF{MU} end

Base.:+(nodof::AbstractNoDOF) = nodof
Base.:+(::IntegerLike, nodof::AbstractNoDOF) = nodof
Base.:+(nodof::AbstractNoDOF, ::IntegerLike) = nodof
Base.:+(nodof::AbstractNoDOF, ::AbstractNoDOF) = nodof

Base.:*(nodof::AbstractNoDOF) = nodof
Base.:*(::IntegerLike, nodof::AbstractNoDOF) = nodof
Base.:*(nodof::AbstractNoDOF, ::IntegerLike) = nodof
Base.:*(nodof::AbstractNoDOF, ::AbstractNoDOF) = nodof


"""
    MeasureBase.NoDOF{MU} <: AbstractNoDOF{MU}

Indicates that there is no way to compute degrees of freedom of a measure
of type `MU` with the given information, e.g. because the DOF are not
a global property of the measure.
"""
struct NoDOF{MU} <: AbstractNoDOF{MU} end


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
export getdof

# Prevent infinite recursion:
@inline _default_getdof(::Type{MU}, ::MU) where {MU} = NoDOF{MU}()
@inline _default_getdof(::Type{MU}, mu_base) where {MU} = getdof(mu_base)

@inline getdof(μ::MU) where {MU} = _default_getdof(MU, basemeasure(μ))


"""
    MeasureBase.NoFastDOF{MU} <: AbstractNoDOF{MU}

Indicates that there is no way to compute degrees of freedom of a measure
of type `MU` with the given information, e.g. because the DOF are not
a global property of the measure.
"""
struct NoFastDOF{MU} <: AbstractNoDOF{MU} end


"""
    fast_dof(μ::MU)

Returns the effective number of degrees of freedom of variates of
measure `μ`, if it can be computed efficiently, otherwise
returns [`NoFastDOF{MU}()`](@ref).

Defaults to `getdof(μ)` and should be specialized for measures for
wich DOF can't be computed instantly.

The effective NDOF my differ from the length of the variates. For example,
the effective NDOF for a Dirichlet distribution with variates of length `n`
is `n - 1`.

Also see [`check_dof`](@ref).
"""
function fast_dof end
export fast_dof

fast_dof(μ) = getdof(μ)


"""
    MeasureBase.some_dof(μ::AbstractMeasure)

Get the DOF at some unspecified point of measure `μ`.

Use with caution!

In general, use [`getdof(μ)`](@ref) instead. `some_dof` is useful
for measures are expected to have a constant DOF of their whole
space but for which there is no way to compute it (or prove that
the DOF is constant of the measurable space).
"""
function some_dof end

function some_dof(μ)
    m = asmeasure(μ)
    _try_direct_dof(m, getdof(m))
end

_try_direct_dof(::AbstractMeasure, dof::IntegerLike) = dof
_try_direct_dof(μ::AbstractMeasure, ::AbstractNoDOF) = _try_local_dof(μ::AbstractMeasure, some_dof(_some_localmeasure(μ)))
_try_local_dof(::AbstractMeasure, dof::IntegerLike) = dof
_try_local_dof(μ::AbstractMeasure, ::AbstractNoDOF) = throw(ArgumentError("Can't determine DOF for measure of type $(nameof(typeof(μ)))"))

_some_localmeasure(μ::AbstractMeasure) = localmeasure(μ, testvalue(μ))


"""
    MeasureBase.check_dof(ν, μ)::Nothing

Check if `ν` and `μ` have the same effective number of degrees of freedom
according to [`MeasureBase.fast_dof`](@ref).
"""
function check_dof end

function check_dof(ν, μ)
    n_ν = fast_dof(ν)
    n_μ = fast_dof(μ)
    # TODO: How to handle this better if DOF is unclear e.g. for HierarchicalMeasures?
    if n_ν isa AbstractNoDOF || n_μ isa AbstractNoDOF
        return true
    end
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
