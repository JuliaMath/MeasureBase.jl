# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# The distribution types `Uniform` and `Normal` stand for the standard
# measures `StdUniform` and `StdNormal` as transport partners:
const _AnyStdDistribution = Union{Uniform,Normal}

_stdmeasure_type(::Type{<:Uniform}) = StdUniform
_stdmeasure_type(::Type{<:Normal}) = StdNormal

MeasureBase.transport_to(::Type{NU}, μ) where {NU<:_AnyStdDistribution} = transport_to(_stdmeasure_type(NU), μ)
MeasureBase.transport_to(ν, ::Type{MU}) where {MU<:_AnyStdDistribution} = transport_to(ν, _stdmeasure_type(MU))

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
