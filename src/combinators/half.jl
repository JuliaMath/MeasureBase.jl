export Half

"""
    Half(μ)

The distribution of `abs(x)` for `x` drawn from `μ`, for measures `μ` with
scalar variates that are symmetric about zero and have no atoms (the fold
doubles densities, which would double the mass of an atom at zero).
"""
struct Half{M} <: AbstractMeasure
    parent::M

    function Half(parent::M) where {M}
        _check_half_parent(parent, mspace_ndims(M), _symmetric_about_zero(parent))
        _atom_free(rootmeasure(parent)) || _throw_half_parent(parent, "a root measure without atoms")
        return new{M}(parent)
    end
end

@inline _check_half_parent(::Any, ::StaticInteger{0}, ::Any) = nothing
@inline _check_half_parent(μ, n::Integer, s) = n == 0 ? _check_half_parent(μ, static(0), s) : _throw_half_parent(μ, "scalar variates")
@inline _check_half_parent(μ, ::NoMSpaceElementSize, ::Any) = _throw_half_parent(μ, "a declared variate rank of zero")
@inline _check_half_parent(μ, ::StaticInteger{0}, ::False) = _throw_half_parent(μ, "symmetry about zero")
@noinline _throw_half_parent(μ, what) =
    throw(ArgumentError("Half requires measures with $what, but got a measure of type $(nameof(typeof(μ)))"))

# Whether a root measure has no atoms, true unless known otherwise:
_atom_free(::Any) = true
_atom_free(::Union{Dirac,Counting,CountingBase,SpikeMixture}) = false
_atom_free(μ::SuperpositionMeasure) = all(_atom_free, μ.components)

# Whether a measure is symmetric about zero, `nothing` if not known:
_symmetric_about_zero(::Any) = nothing
_symmetric_about_zero(::Union{StdNormal,StdLogistic}) = static(true)
_symmetric_about_zero(::Union{StdExponential,StdUniform}) = static(false)
_symmetric_about_zero(μ::WeightedMeasure) = _symmetric_about_zero(μ.base)

@inline mspace_elsize(μ::Half) = mspace_elsize(μ.parent)
@inline mspace_flatsize(μ::Half) = mspace_flatsize(μ.parent)
@inline mspace_flatsize(::Type{<:Half{M}}) where {M} = mspace_flatsize(M)
@inline mspace_ndims(::Type{<:Half{M}}) where {M} = mspace_ndims(M)
@inline stream_ndims(::Type{<:Half{M}}) where {M} = stream_ndims(M)
@inline fixed_stream_size(::Type{<:Half{M}}) where {M} = fixed_stream_size(M)
@inline preferred_stdmeasure(::Type{<:Half}) = StdUniform

function Base.show(io::IO, μ::Half)
    print(io, "Half")
    show(io, μ.parent)
end

unhalf(μ::Half) = μ.parent

@inline function basemeasure(μ::Half)
    weightedmeasure(logtwo, basemeasure(unhalf(μ)))
end

@inline rand_impl(ctx::GenContext, μ::Half) = abs(rand_impl(ctx, unhalf(μ)))
@inline batched_rand_impl(ctx::GenContext, μ::Half, sz::SizeLike) = abs.(batched_rand_impl(ctx, unhalf(μ), sz))

function logdensityof_impl(μ::Half, x)
    ld = logdensityof(unhalf(μ), x) + logtwo
    return ifelse(x ≥ 0, ld, oftype(ld, -Inf))
end

logdensity_def(μ::Half, x) = logdensity_def(unhalf(μ), x)

@inline insupportof(d::Half, x) = _insupport_and(x ≥ 0, insupportof(unhalf(d), x))

testvalue(::Type{T}, ::Half) where {T} = one(T)

massof(μ::Half) = massof(unhalf(μ))

# The fold doubles the upper half of the parent, the Stieltjes measure
# functions scale with the mass of the parent (`smf(μ.parent, 0)` is half
# of it):
smf(μ::Half, x) = _half_smf(μ.parent, x)
csmf(μ::Half, x) = 2 * csmf(μ.parent, max(x, zero(x)))
invsmf(μ::Half, p) = _nan_outside_mass(μ, p, invsmf(μ.parent, _half_smf_arg(μ.parent, p)))
invcsmf(μ::Half, p) = _nan_outside_mass(μ, p, invcsmf(μ.parent, p / 2))

@inline _half_smf(ν, x) = 2 * (smf(ν, max(x, zero(x))) - smf(ν, zero(x)))
@inline _half_smf_arg(ν, p) = smf(ν, zero(p)) + p / 2

@inline _nan_outside_mass(μ, p, x) = ifelse((zero(p) <= p) & (p <= massof(μ)), x, oftype(x, NaN))

# Transports map the normalized law, the fold of the normalized parent:
@inline transport_to_std(::Type{StdUniform}, μ::Half, x) = _nan_outside(μ, x, _half_smf(normalize(μ.parent), x))
@inline function transport_from_std(::Type{StdUniform}, μ::Half, u)
    ν = normalize(μ.parent)
    _nan_outside(StdUniform(), u, invsmf(ν, _unit_interior(_half_smf_arg(ν, u))))
end
