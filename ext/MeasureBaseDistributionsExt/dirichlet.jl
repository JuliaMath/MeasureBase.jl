# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

const DirichletMeasure = AsMeasure{<:Dirichlet}

MeasureBase.getdof(d::Dirichlet) = length(d) - 1
MeasureBase.getdof(m::DirichletMeasure) = getdof(m.obj)

MeasureBase.transport_origin(d::Dirichlet) = StdUniform()^getdof(d)



function _dirichlet_beta_trafo(α::Real, β::Real, x::Real)
    R = float(promote_type(typeof(α), typeof(β), typeof(x)))
    convert(R, transport_def(Beta(α, β), StdUniform(), x))::R
end

_a_times_one_minus_b(a::Real, b::Real) = a * (1 - b)

function MeasureBase.from_origin(ν::Dirichlet, x)
    # See M. J. Betancourt, "Cruising The Simplex: Hamiltonian Monte Carlo and the Dirichlet Distribution",
    # https://arxiv.org/abs/1010.3436

    @_adignore @argcheck length(ν) == length(x) + 1

    αs = _dropfront(_rev_cumsum(ν.alpha))
    βs = _dropback(ν.alpha)
    beta_v = _fwddiff(_dirichlet_beta_trafo).(αs, βs, x)
    beta_v_cp = _exp_cumsum_log(_pushfront(beta_v, 1))
    beta_v_ext = _pushback(beta_v, 0)
    _fwddiff(_a_times_one_minus_b).(beta_v_cp, beta_v_ext)
end


function _inv_dirichlet_beta_trafo(α::Real, β::Real, beta_v::Real)
    R = float(promote_type(typeof(α), typeof(β), typeof(beta_v)))
    convert(R, transport_def(StdUniform(), Beta(α, β), beta_v))::R
end

# ToDo: Find efficient pullback for this:
function _dirichlet_variate_to_beta_v(y::AbstractVector{<:Real})
    beta_v = similar(y, length(eachindex(y)) - 1)
    @assert firstindex(beta_v) == firstindex(y)
    @assert lastindex(beta_v) == lastindex(y) - 1
    T = eltype(y)
    sum_log_beta_v::T = 0
    @inbounds for i in eachindex(beta_v)
        beta_v[i] = 1 - y[i] / exp(sum_log_beta_v)
        sum_log_beta_v += log(beta_v[i])
    end
    return beta_v
end

function MeasureBase.to_origin(ν::Dirichlet, y)
    @_adignore @argcheck length(ν) == length(y)
    αs = _dropfront(_rev_cumsum(ν.alpha))
    βs = _dropback(ν.alpha)
    beta_v = _dirichlet_variate_to_beta_v(y)
    _fwddiff(_inv_dirichlet_beta_trafo).(αs, βs, beta_v)
end
