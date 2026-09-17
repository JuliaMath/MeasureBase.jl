# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

const _AnyStdUniform = Union{StandardUniform,Uniform}
const _AnyStdNormal = Union{StandardNormal,Normal}

const _AnyStdDistribution = Union{_AnyStdUniform,_AnyStdNormal}

_std_dist(::Type{<:_AnyStdUniform}) = StandardUniform
_std_dist(::Type{<:_AnyStdNormal}) = StandardNormal

_std_dist(::Type{D}, ::StaticInt{1}) where {D<:_AnyStdDistribution} = D()
_std_dist(::Type{D}, dof) where {D<:_AnyStdDistribution} = D(dynamic(dof))
_std_dist_for(::Type{D}, μ::Any) where {D<:_AnyStdDistribution} = _std_dist(_std_dist(D), getdof(μ))

MeasureBase.transport_to(::Type{NU}, μ) where {NU<:_AnyStdDistribution} = transport_to(_std_dist_for(NU, μ), μ)
MeasureBase.transport_to(ν, ::Type{MU}) where {MU<:_AnyStdDistribution} = transport_to(ν, _std_dist_for(MU, ν))

# Disambiguation between the type forms of standard measures and distributions:
function MeasureBase.transport_to(::Type{NU}, ::Type{MU}) where {NU<:_AnyStdDistribution,MU<:_AnyStdDistribution}
    _throw_two_std_types()
end
function MeasureBase.transport_to(::Type{NU}, ::Type{MU}) where {NU<:StdMeasure,MU<:_AnyStdDistribution}
    _throw_two_std_types()
end
function MeasureBase.transport_to(::Type{NU}, ::Type{MU}) where {NU<:_AnyStdDistribution,MU<:StdMeasure}
    _throw_two_std_types()
end
function _throw_two_std_types()
    throw(
        ArgumentError(
            "Can't construct a transport function between the types of two standard measures, need a measure instance on one side",
        ),
    )
end
