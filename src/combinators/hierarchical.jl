export HierarchicalMeasure


# TODO: Document and use FlattenMode
abstract type FlattenMode end
struct NoFlatten <: FlattenMode end
struct AutoFlatten <: FlattenMode end


struct HierarchicalMeasure{F,M<:AbstractMeasure,FM<:FlattenMode} <: AbstractMeasure
    f::F
    m::M
    flatten_mode::FM
end

# TODO: Document
const HierarchicalProductMeasure{F,M<:AbstractMeasure} = HierarchicalMeasure{F,M,NoFlatten}
export HierarchicalProductMeasure

HierarchicalProductMeasure(f, m::AbstractMeasure) = HierarchicalMeasure(f, m, NoFlatten())

# TODO: Document
const FlatHierarchicalMeasure{F,M<:AbstractMeasure} = HierarchicalMeasure{F,M,AutoFlatten}
export FlatHierarchicalMeasure

FlatHierarchicalMeasure(f, m::AbstractMeasure) = HierarchicalMeasure(f, m, AutoFlatten())

HierarchicalMeasure(f, m::AbstractMeasure) = FlatHierarchicalMeasure(f, m)



function _split_variate_after(::NoFlatten, μ::AbstractMeasure, x::Tuple{2})
    @assert x isa Tuple{2}
    return x[1], x[2]
end


function _split_variate_after(::AutoFlatten, μ::AbstractMeasure, x)
    a_test = testvalue(μ)
    return _autosplit_variate_after_testvalue(a_test, x)
end

function _autosplit_variate_after_testvalue(::Any, x)
    @assert x isa Tuple{2}
    return x[1], x[2]
end

function _autosplit_variate_after_testvalue(a_test::AbstractVector, x::AbstractVector)
    n, m = length(eachindex(a_test)), length(eachindex(x))
    # TODO: Use getindex or view?
    return x[begin:n], x[begin+n:m]
end

function _autosplit_variate_after_testvalue(::Tuple{N}, x::Tuple{M}) where {N,M}
    return ntuple(i -> x[i], Val(1:N)), ntuple(i -> x[i], Val(N+1:M))
end

@generated function _autosplit_variate_after_testvalue(::NamedTuple{names_a}, x::NamedTuple{names}) where {names_a,names}
    # TODO: implement
    @assert false
end



_combine_variates(::NoFlatten, a::Any, b::Any) = (a, b)


_combine_variates(::AutoFlatten, a::Any, b::Any) = _autoflat_combine_variates(a, b)

_autoflat_combine_variates(a::Any, b::Any) = (a, b)

_autoflat_combine_variates(a::AbstractVector, b::AbstractVector) = vcat(a, b)

_autoflat_combine_variates(a::Tuple, b::Tuple) = (a, b)

# TODO: Check that names don't overlap:
_autoflat_combine_variates(a::NamedTuple, b::NamedTuple) = merge(a, b)


_local_productmeasure(::NoFlatten, μ1, μ2) = productmeasure(μ1, μ2)

# TODO: _local_productmeasure(::AutoFlatten, μ1, μ2) = productmeasure(μ1, μ2)
# Needs a FlatProductMeasure type.

function _localmeasure_with_rest(μ::HierarchicalProductMeasure, x)
    μ_primary = μ.m
    local_primary, x_secondary = _localmeasure_with_rest(μ_primary, x)
    μ_secondary = μ.f(x_secondary)
    local_secondary, x_rest = _localmeasure_with_rest(μ_secondary, x_secondary)
    return _local_productmeasure(μ.flatten_mode, local_primary, local_secondary), x_rest
end

function _localmeasure_with_rest(μ::AbstractMeasure, x)
    x_checked = checked_arg(μ, x)
    return localmeasure(μ, x_checked), Fill(zero(eltype(x)), 0)
end

function localmeasure(μ::HierarchicalProductMeasure, x)
    h_local, x_rest = _localmeasure_with_rest(μ, x)
    if !isempty(x_rest)
        throw(ArgumentError("Variate too long while computing localmeasure of HierarchicalMeasure"))
    end
    return h_local
end


@inline insupport(::HierarchicalMeasure, x) = NoFastInsupport()

@inline getdof(μ::HierarchicalMeasure) = NoDOF{typeof(μ)}()

# Bypass `checked_arg`, would require potentially costly evaluation of h.f:
@inline checked_arg(::HierarchicalMeasure, x) = x

rootmeasure(::HierarchicalMeasure) = throw(ArgumentError("root measure is implicit, but can't be instantiated, for HierarchicalMeasure"))

basemeasure(::HierarchicalMeasure) = throw(ArgumentError("basemeasure is not available for HierarchicalMeasure"))

logdensity_def(::HierarchicalMeasure, x) = throw(ArgumentError("logdensity_def is not available for HierarchicalMeasure"))


# # TODO: Default implementation of unsafe_logdensityof is a bit inefficient
# # for AutoFlatten, since variate will be split in `localmeasure` and then
# # split again in log-density evaluation. Maybe add something like
# function unsafe_logdensityof(h::HierarchicalMeasure, x)
#     local_primary, local_secondary, x_primary, x_secondary = ...
#     # Need to call full logdensityof for h_secondary since x_secondary hasn't
#     # been checked yet:
#     unsafe_logdensityof(local_primary, x_primary) + logdensityof(local_secondary, x_secondary)
# end


function Base.rand(rng::Random.AbstractRNG, ::Type{T}, h::HierarchicalMeasure) where {T<:Real}
    x_primary = rand(rng, T, h.m)
    x_secondary = rand(rng, T, h.f(x_primary))
    return _combine_variates(h.flatten_mode, x_primary, x_secondary)
end



function _to_std_with_rest(flatten_mode::FlattenMode, ν_inner::StdMeasure, μ::HierarchicalMeasure, x)
    μ_primary = μ.m
    y_primary, x_secondary = _to_std_with_rest(flatten_mode, ν_inner, μ_primary, x)
    μ_secondary = μ.f(x_secondary)
    y_secondary, x_rest = _to_std_with_rest(flatten_mode, ν_inner, μ_secondary, x_secondary)
    return _combine_variates(μ.flatten_mode, y_primary, y_secondary), x_rest
end

function _to_std_with_rest(flatten_mode::FlattenMode, ν_inner::StdMeasure, μ::AbstractMeasure, x)
    dof_μ = getdof(μ)
    x_μ, x_rest = _split_variate_after(flatten_mode, μ, x)
    y = transport_to(ν_inner^dof_μ, μ, x_μ)
    return y, x_rest
end

function transport_def(ν::_PowerStdMeasure{1}, μ::HierarchicalMeasure, x)
    ν_inner = _get_inner_stdmeasure(ν)
    y, x_rest = _to_std_with_rest(ν_inner, μ, x)
    if !isempty(x_rest)
        throw(ArgumentError("Variate too long during transport involving HierarchicalMeasure"))
    end
    return y
end


function _from_std_with_rest(ν::HierarchicalMeasure, μ_inner::StdMeasure, x)
    ν_primary = ν.m
    y_primary, x_secondary = _from_std_with_rest(ν_primary, μ_inner, x)
    ν_secondary = ν.f(y_primary)
    y_secondary, x_rest = _from_std_with_rest(ν_secondary, μ_inner, x_secondary)
    return _combine_variates(ν.flatten_mode, y_primary, y_secondary), x_rest
end

function _from_std_with_rest(ν::AbstractMeasure, μ_inner::StdMeasure, x)
    dof_ν = getdof(ν)
    len_x = length(eachindex(x))

    # Since we can't check DOF of original HierarchicalMeasure, we could "run out x" if
    # the original x was too short. `transport_to` below will detect this, but better
    # throw a more informative exception here:
    if len_x < dof_ν
        throw(ArgumentError("Variate too short during transport involving HierarchicalMeasure"))
    end

    y = transport_to(ν, μ_inner^dof_ν, x[begin:begin+dof_ν-1])
    x_rest = Fill(zero(eltype(x)), dof_ν - len_x)
    return y, x_rest
end

function transport_def(ν::HierarchicalMeasure, μ::_PowerStdMeasure{1}, x)
    # Sanity check, should be checked by transport machinery already:
    @assert getdof(μ) == length(eachindex(x)) && x isa AbstractVector
    μ_inner = _get_inner_stdmeasure(μ)
    y, x_rest = _from_std_with_rest(ν, μ_inner, x)
    if !isempty(x_rest)
        throw(ArgumentError("Variate too long during transport involving HierarchicalMeasure"))
    end
    return y
end
