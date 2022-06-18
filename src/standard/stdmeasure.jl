abstract type StdMeasure <: AbstractMeasure end

StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeof(randexp)) = StdExponential()


@inline check_dof(::StdMeasure, ::StdMeasure) = nothing


@inline vartransform_def(::MU, μ::MU, x) where {MU<:StdMeasure} = x

function vartransform_def(ν::StdMeasure, μ::PowerMeasure{<:StdMeasure}, x)
    return vartransform_def(ν, μ.parent, only(x))
end

function vartransform_def(ν::PowerMeasure{<:StdMeasure}, μ::StdMeasure, x)
    return Fill(vartransform_def(ν.parent, μ, only(x)), map(length, ν.axes)...)
end

function vartransform_def(ν::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}}, μ::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}}, x)
    return vartransform(ν.parent, μ.parent).(x)
end

function vartransform_def(ν::PowerMeasure{<:StdMeasure,<:NTuple{N,Base.OneTo}}, μ::PowerMeasure{<:StdMeasure,<:NTuple{M,Base.OneTo}}, x) where {N,M}
    return reshape(vartransform(ν.parent, μ.parent).(x), map(length, ν.axes)...)
end


# Implement vartransform(NU::Type{<:StdMeasure}, μ) and vartransform(ν, MU::Type{<:StdMeasure}):

_std_measure(::Type{M}, ::StaticInt{1}) where {M<:StdMeasure} = M()
_std_measure(::Type{M}, dof::Integer) where {M<:StdMeasure} = M()^dof
_std_measure_for(::Type{M}, μ::Any) where {M<:StdMeasure} = _std_measure(M, getdof(μ))

MeasureBase.vartransform(::Type{NU}, μ) where {NU<:StdMeasure} = vartransform(_std_measure_for(NU, μ), μ)
MeasureBase.vartransform(ν, ::Type{MU}) where {MU<:StdMeasure} = vartransform(ν, _std_measure_for(MU, ν))
