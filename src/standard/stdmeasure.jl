abstract type StdMeasure <: AbstractMeasure end

StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeof(randexp)) = StdExponential()
StdMeasure(::typeof(randn)) = StdNormal()

"""
    MeasureBase.StdPowerMeasure{MU<:StdMeasure,N}

The type of an `N`-dimensional power of a standard measure of type `MU`.
"""
const StdPowerMeasure{MU<:StdMeasure,N} = PowerMeasure{MU,<:NTuple{N,OneToLike}}

@inline mspace_elsize(::StdMeasure) = ()

@inline check_dof(::StdMeasure, ::StdMeasure) = nothing

@inline transport_def(::MU, μ::MU, x) where {MU<:StdMeasure} = x
