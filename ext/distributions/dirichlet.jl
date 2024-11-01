# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

const DirichletMeasure = AsMeasure{<:Dirichlet}

MeasureBase.getdof(m::DirichletMeasure) = length(m.obj) - 1

MeasureBase.transport_origin(m::DirichletMeasure) = StdUniform()^getdof(m)



function _dirichlet_beta_trafo(α::Real, β::Real, x::Real)
    R = float(promote_type(typeof(α), typeof(β), typeof(x)))
    convert(R, transport_def(Beta(α, β), StdUniform(), x))::R
end

_a_times_one_minus_b(a::Real, b::Real) = a * (1 - b)

function MeasureBase.from_origin(ν::Dirichlet, x)
    # See M. J. Betancourt, "Cruising The Simplex: Hamiltonian Monte Carlo and the Dirichlet Distribution",
    # https://arxiv.org/abs/1010.3436

    # Sanity check (TODO - remove?):
    @_adignore @argcheck length(ν) == length(x) + 1

    αs = _dropfront(_rev_cumsum(ν.alpha))
    βs = _dropback(ν.alpha)
    beta_v = fwddiff(_dirichlet_beta_trafo).(αs, βs, x)
    beta_v_cp = _exp_cumsum_log(_pushfront(beta_v, 1))
    beta_v_ext = _pushback(beta_v, 0)
    fwddiff(_a_times_one_minus_b).(beta_v_cp, beta_v_ext)
end

# ToDo: MeasureBase.to_origin(ν::Dirichlet, y)
