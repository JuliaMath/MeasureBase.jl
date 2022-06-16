abstract type StdMeasure <: AbstractMeasure end

StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeof(randexp)) = StdExponential()


@inline vartransform_def(::MU, ::MU, x::Real) where {MU<:StdMeasure} = x

function vartransform_def(ν::StdMeasure, μ::PowerMeasure{<:StdMeasure}, x)
    check_dof(ν, μ)
    vartransform_def(ν, μ.parent, only(x))
end

function vartransform_def(ν::PowerMeasure{<:StdMeasure}, μ::StdMeasure, x)
    check_dof(ν, μ)
    Fill(vartransform_def(ν.parent, μ, only(x)), map(length, ν.axes)...)
end

function vartransform_def(ν::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}}, μ::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}}, x)
    check_dof(ν, μ)
    vartransform(ν.parent, μ.parent).(x)
end

function vartransform_def(ν::PowerMeasure{<:StdMeasure,<:NTuple{N,Base.OneTo}}, μ::PowerMeasure{<:StdMeasure,<:NTuple{M,Base.OneTo}}, x) where {N,M}
    check_dof(ν, μ)
    reshape(vartransform(ν.parent, μ.parent).(x), map(length, ν.axes)...)
end
