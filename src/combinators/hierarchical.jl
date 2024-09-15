export HierarchicalMeasure


struct HierarchicalMeasure{F,M<:AbstractMeasure} <: AbstractMeasure
    f::F
    m::M
    dof_m::Int
end


function HierarchicalMeasure(f, m::AbstractMeasure, ::NoDOF)
    throw(ArgumentError("Primary measure in HierarchicalMeasure must have fixed and known DOF"))
end

HierarchicalMeasure(f, m::AbstractMeasure) = HierarchicalMeasure(f, m, dynamic(getdof(m)))


function _split_variate(h::HierarchicalMeasure, x)
    # TODO: Splitting x will be more complicated in general:
    x_primary, x_secondary = x
    return (x_primary, x_secondary)
end


function _combine_variates(x_primary, x_secondary)
    # TODO: Must offer optional flattening
    return (x_primary, x_secondary)
end


function local_measure(h::HierarchicalMeasure, x)
    x_primary, x_secondary = _split_variate(h, x)
    m_primary = h.m
    m_primary_local = local_measure(m_primary, x_primary)
    m_secondary = m.f(x_secondary)
    m_secondary_local = local_measure(m_secondary, x_secondary)
    # TODO: Must optionally return a flattened product measure
    return productmeasure(m_primary_local, m_secondary_local)
end


@inline function insupport(h::HierarchicalMeasure, x)
    # Only test primary for efficiency:
    x_primary = _split_variate(h, x)[1]
    insupport(h.m, x_primary)
end


#!!!!!!! WON'T WORK: Only use primary measure for efficiency:
logdensity_type(h::HierarchicalMeasure{F,M}, ::Type{T}) where {F,M,T} = unstatic(float(logdensity_type(M, T)))

# Can't implement logdensity_def(::HierarchicalMeasure, x) directly.

# Can't implement getdof(::HierarchicalMeasure) efficiently

# No way to return a functional base measure:
struct _BaseOfHierarchicalMeasure{F,M<:AbstractMeasure} <: AbstractMeasure end
@inline basemeasure(::HierarchicalMeasure{F,M}) where {F,M} = _BaseOfHierarchicalMeasure{F,M}()

@inline getdof(μ::HierarchicalMeasure) = NoDOF{typeof(μ)}()

# Bypass `checked_arg`, would require potentially costly evaluation of h.f:
@inline checked_arg(::HierarchicalMeasure, x) = x

function unsafe_logdensityof(h::HierarchicalMeasure, x)
    x_primary, x_secondary = _split_variate(h, x)
    h_primary, h_secondary = h.m, h.f(x_secondary)
    unsafe_logdensityof(h_primary, x_primary) + logdensityof(h_secondary, x_secondary)
end


function Base.rand(rng::Random.AbstractRNG, ::Type{T}, h::HierarchicalMeasure) where {T<:Real}
    x_primary = rand(rng, T, h.m)
    x_secondary = rand(rng, T, h.f(x_primary))
    return _combine_variates(x_primary, x_secondary)
end


function _split_measure_at(μ::PowerMeasure{M, Tuple{R}}, n::Integer) where {M<:StdMeasure,R}
    dof_μ = getdof(μ)
    return M()^n, M()^(dof_μ - n)
end


function transport_def(
    ν::PowerMeasure{M, Tuple{R}},
    μ::HierarchicalMeasure,
    x,
) where {M<:StdMeasure,R}
    ν_primary, ν_secondary = _split_measure_at(ν, μ.dof_m)
    x_primary, x_secondary = _split_variate(μ, x)
    μ_primary = μ.m
    μ_secondary = μ.f(x_secondary)
    y_primary = transport_to(ν_primary, μ_primary, x_primary)
    y_secondary = transport_to(ν_secondary, μ_secondary, x_secondary)
    return vcat(y_primary, y_secondary)
end


function transport_def(
    ν::HierarchicalMeasure,
    μ::PowerMeasure{M, Tuple{R}},
    x,
) where {M<:StdMeasure,R}
    dof_primary = ν.dof_m
    μ_primary, μ_secondary = _split_measure_at(μ, dof_primary)
    x_primary, x_secondary = x[begin:begin+dof_primary-1], x[begin+dof_primary:end]
    ν_primary = ν.m
    y_primary = transport_to(ν_primary, μ_primary, x_primary)
    ν_secondary = ν.f(y_primary)
    y_secondary = transport_to(ν_secondary, μ_secondary, x_secondary)
    return _combine_variates(y_primary, y_secondary)
end
