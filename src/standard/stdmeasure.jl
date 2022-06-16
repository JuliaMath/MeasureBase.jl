abstract type StdMeasure <: AbstractMeasure end

StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeof(randexp)) = StdExponential()


@inline check_dof(::StdMeasure, ::StdMeasure) = nothing

@inline checked_var(::StdMeasure, x::Real) = x

@propagate_inbounds function checked_var(::StdMeasure, x::Any)
    @boundscheck throw(ArgumentError("Invalid variate type for measure"))
end


@inline vartransform_def(::MU, μ::MU, x) where {MU<:StdMeasure} = x


function vartransform_def(ν::StdMeasure, μ::PowerMeasure{<:StdMeasure}, x)
    check_dof(ν, μ)
    return vartransform_def(ν, μ.parent, only(x))
end

function vartransform_def(ν::PowerMeasure{<:StdMeasure}, μ::StdMeasure, x)
    check_dof(ν, μ)
    return Fill(vartransform_def(ν.parent, μ, only(x)), map(length, ν.axes)...)
end

function vartransform_def(ν::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}}, μ::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}}, x)
    check_dof(ν, μ)
    return vartransform(ν.parent, μ.parent).(x)
end

function vartransform_def(ν::PowerMeasure{<:StdMeasure,<:NTuple{N,Base.OneTo}}, μ::PowerMeasure{<:StdMeasure,<:NTuple{M,Base.OneTo}}, x) where {N,M}
    check_dof(ν, μ)
    return reshape(vartransform(ν.parent, μ.parent).(x), map(length, ν.axes)...)
end
