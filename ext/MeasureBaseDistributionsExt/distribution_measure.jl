# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).


const DistributionMeasure{F<:VariateForm,S<:ValueSupport,D<:Distribution{F,S}} = AsMeasure{D}

@inline MeasureBase.AbstractMeasure(obj::Distribution) = AsMeasure{typeof(obj)}(obj)
@inline Base.convert(::Type{AbstractMeasure}, obj::Distribution) = AbstractMeasure(obj)

@inline Distributions.Distribution(m::DistributionMeasure) = m.obj
@inline Distributions.Distribution{F}(m::DistributionMeasure{F}) where {F<:VariateForm} = Distribution(m)
@inline Distributions.Distribution{F,S}(m::DistributionMeasure{F,S}) where {F<:VariateForm,S<:ValueSupport} = Distribution(m)

@inline Base.convert(::Type{Distribution}, m::DistributionMeasure) = Distribution(m)
@inline Base.convert(::Type{Distribution{F}}, m::DistributionMeasure{F}) where {F<:VariateForm} = Distribution(m)
@inline Base.convert(::Type{Distribution{F,S}}, m::DistributionMeasure{F,S}) where {F<:VariateForm,S<:ValueSupport} = Distribution(m)


# Distributions' samplers run on host compute units and draw in the
# precision of the distribution parameters. Where that is lower than the
# precision of the context, and on other compute units, variates are
# generated from standard variates of the context precision via the
# transports. Draws of continuous distributions without standard transports
# are converted to the context precision afterwards, only their type
# follows the context then. Atoms of discrete distributions stay exact.
# On other compute units, distributions without standard transports are
# sampled on the host, see `_host_rand`.
function MeasureBase.rand_impl(ctx::GenContext, m::DistributionMeasure)
    d = m.obj
    if _native_rand(ctx, d)
        _native_draw(get_precision(ctx), d, _dist_rand(get_rng(ctx), d))
    else
        _nonnative_rand(ctx, m, ())
    end
end

function MeasureBase.batched_rand_impl(ctx::GenContext, m::DistributionMeasure, sz::SizeLike)
    d = m.obj
    if _native_rand(ctx, d)
        _flat_powrand(get_rng(ctx), get_precision(ctx), d, sz)
    else
        _nonnative_rand(ctx, m, sz)
    end
end

function _nonnative_rand(ctx::GenContext, m::DistributionMeasure, sz::SizeLike)
    if _has_std_transport(m.obj)
        MeasureBase._rand_default(ctx, m, sz, MeasureBase._NoRandImpl())
    else
        _host_rand(ctx, m, sz)
    end
end

# Distributions' samplers draw array-like variates on the host, in a host
# context derived from the context (whose RNG advances), and the variates
# are moved to the compute unit:
function _host_rand(ctx::GenContext, m::DistributionMeasure{<:ArrayLikeVariate}, sz::SizeLike)
    MeasureBase._can_derive_host_context(ctx) || _throw_host_only_sampler(m.obj, "it can't run within compiled code")
    host_ctx = MeasureBase._host_context(ctx)
    MeasureBase._on_compute_unit(ctx, MeasureBase.batched_rand_impl(host_ctx, m, asnonstatic(sz)))
end
_host_rand(::GenContext, m::DistributionMeasure, ::SizeLike) =
    _throw_host_only_sampler(m.obj, "its variates can't be moved to other compute units")
# Bernoulli variates compare uniform variates with the success probability,
# branch-free and exact:
function _nonnative_rand(ctx::GenContext, m::AsMeasure{<:Bernoulli}, sz::SizeLike)
    MeasureBase._rand_bulk(ctx, sz) .< convert(get_precision(ctx), Distributions.succprob(m.obj))
end

@inline function _native_rand(ctx::GenContext, d::Distribution)
    MeasureBase._on_host(ctx) && (_native_precision_suffices(get_precision(ctx), d) || !_has_std_transport(d))
end
@inline function _native_precision_suffices(::Type{T}, d::Distribution) where {T<:Real}
    P = float(Distributions.partype(d))
    promote_type(P, T) == P
end
@inline _has_std_transport(d::Distribution) = _has_std_transport(typeof(d))
@inline _has_std_transport(::Type{D}) where {D<:Distribution} = MeasureBase.preferred_stdmeasure(D) <: StdMeasure

# Arrays of distributions without standard transports have no batched
# kernels, so their products draw in host loops (also Bernoulli ones):
@inline MeasureBase._host_drawn_marginals(::Type{AsMeasure{D}}) where {D<:Distribution{<:ArrayLikeVariate}} =
    !_has_std_transport(D)

@noinline _throw_host_only_sampler(d::Distribution, reason::AbstractString) = throw(ArgumentError(
    "$(nameof(typeof(d))) has no standard transport and Distributions' sampler for it runs only on host compute units, $reason"
))

@inline _native_draw(::Type{T}, ::Distribution{<:ArrayLikeVariate,Continuous}, x) where {T<:Real} = convert_realtype(T, x)
@inline _native_draw(::Type{<:Real}, ::Distribution, x) = x

# Transformed atoms stay exact as well:
@inline MeasureBase._pushfwd_draw(::GenContext, ::DistributionMeasure{<:Any,Discrete}, y) = y

# Constant RNGs (for test values) only draw floats, but Distributions'
# samplers draw integer ranges and may reject constant draws forever, so
# univariate distributions are sampled by inversion then (in the type of
# their native draws), multinomials by the inversion of their sequential
# binomial draws:
@inline _dist_rand(rng::AbstractRNG, d) = rand(rng, d)
@inline _dist_rand(rng::ConstantRNG, d::Distribution{Univariate}) = convert(eltype(d), Distributions.quantile(d, rand(rng)))
# Samplers of non-array variates (e.g. LKJCholesky) fail with constant draws,
# they draw a fixed variate:
@inline _dist_rand(::ConstantRNG, d::Distribution{<:Union{Distributions.CholeskyVariate,Distributions.NamedTupleVariate}}) =
    rand(Xoshiro(0), d)
# Those of correlation matrices and of the sphere give variates outside of
# the support with constant draws, they draw a fixed variate as well:
@inline _dist_rand(::ConstantRNG, d::Union{Distributions.LKJ,Distributions.VonMisesFisher}) = rand(Xoshiro(0), d)
_dist_rand(rng::ConstantRNG, d::Distributions.Multinomial) = _multinomial_by_inversion(rand(rng), Distributions.ntrials(d), Distributions.probs(d))
function _dist_rand(rng::ConstantRNG, d::Distributions.DirichletMultinomial)
    _multinomial_by_inversion(rand(rng), Distributions.ntrials(d), d.α ./ sum(d.α))
end

function _multinomial_by_inversion(u::Real, n::Integer, p::AbstractVector{<:Real})
    x = zeros(Int, length(p))
    n_rem, p_rem = n, one(eltype(p))
    for i in firstindex(p):(lastindex(p) - 1)
        k = Int(Distributions.quantile(Distributions.Binomial(n_rem, clamp(p[i] / p_rem, 0, 1)), u))
        x[i] = k
        n_rem -= k
        p_rem -= p[i]
    end
    x[end] = n_rem
    return x
end

# A single variate for zero batch dimensions, flat batches otherwise, static
# arrays for static sizes of univariate distributions. Constant RNGs give
# constant batches of the single variate:
_flat_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution, sz::SizeLike) where {T<:Real} = _sized_powrand(rng, T, d, sz)
_flat_powrand(rng::ConstantRNG, ::Type{T}, d::Distribution, sz::SizeLike) where {T<:Real} = _const_powrand(T, d, _dist_rand(rng, d), sz)

_sized_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution, ::Tuple{}) where {T<:Real} = _native_draw(T, d, _dist_rand(rng, d))
_sized_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution, sz::SizeLike) where {T<:Real} = _flat_powrand_batch(rng, T, d, asnonstatic(sz))
function _sized_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution{Univariate}, sz::StaticSizeLike) where {T<:Real}
    smp = Distributions.sampler(d)
    A = similar_type(staticarray_type(T, canonical_size(sz)), _native_draw_type(T, d))
    A(ntuple(_ -> _native_draw(T, d, rand(rng, smp)), Val(length(A))))
end
_sized_powrand(rng::AbstractRNG, ::Type{T}, d::Distribution{Univariate}, ::Tuple{}) where {T<:Real} = _native_draw(T, d, _dist_rand(rng, d))

@inline _native_draw_type(::Type{T}, d::Distribution) where {T<:Real} = typeof(_native_draw(T, d, zero(eltype(d))))

_const_powrand(::Type{T}, d::Distribution, x, ::Tuple{}) where {T<:Real} = _native_draw(T, d, x)
_const_powrand(::Type{T}, d::Distribution{Univariate}, x, ::Tuple{}) where {T<:Real} = _native_draw(T, d, x)
_const_powrand(::Type{T}, d::Distribution{Univariate}, x, sz::SizeLike) where {T<:Real} = maybestatic_fill(_native_draw(T, d, x), sz)
function _const_powrand(::Type{T}, d::Distribution, x, sz::SizeLike) where {T<:Real}
    y = _native_draw(T, d, x)
    reshape(repeat(vec(y), prod(asnonstatic(sz))), size(y)..., asnonstatic(sz)...)
end

function _flat_powrand_batch(rng::AbstractRNG, ::Type{T}, d::Distribution{<:ArrayLikeVariate{0}}, sz::Dims) where {T<:Real}
    _native_draw(T, d, reshape(rand(rng, d, prod(sz)), sz...))
end

function _flat_powrand_batch(rng::AbstractRNG, ::Type{T}, d::Distribution{<:ArrayLikeVariate{1}}, sz::Dims) where {T<:Real}
    _native_draw(T, d, reshape(rand(rng, d, prod(sz)), size(d)..., sz...))
end

function _flat_powrand_batch(rng::AbstractRNG, ::Type{T}, d::ReshapedDistribution{N,<:Any,<:Distribution{<:ArrayLikeVariate{1}}}, sz::Dims) where {T<:Real,N}
    _native_draw(T, d, reshape(rand(rng, d.dist, prod(sz)), d.dims..., sz...))
end

function _flat_powrand_batch(rng::AbstractRNG, ::Type{T}, d::Distribution, sz::Dims) where {T<:Real}
    flatview(ArrayOfSimilarArrays(_native_draw(T, d, rand(rng, d, sz))))
end


@inline DensityInterface.densityof(m::DistributionMeasure) = densityof(m.obj)
@inline DensityInterface.logdensityof(m::DistributionMeasure) = logdensityof(m.obj)

@inline MeasureBase.logdensity_def(m::DistributionMeasure, x) = DensityInterface.logdensityof(m.obj, x)

# Distributions' densities may throw outside of the support (e.g.
# InverseGamma below zero, Wishart at matrices that aren't positive
# definite), so they are only evaluated inside of it. Distributions without
# an `insupport` method (Distributions' interface doesn't require one for
# multivariate and matrix-variate distributions) are evaluated at all finite
# variates. Plain conditions branch, traced or device ones select (see
# `MeasureBase._select`). The family kernels don't throw and keep their
# branch-free masks:
@inline function MeasureBase.logdensityof_impl(m::DistributionMeasure, x)
    R = _logd_type(m.obj, x)
    MeasureBase._select(_logpdf_domain(m, x), () -> convert(R, MeasureBase.unsafe_logdensityof(m, x)), () -> R(-Inf))
end
@inline _logd_type(d::Distribution, x) = float(promote_type(_real_eltype(x), Distributions.partype(d)))
@inline _real_eltype(x) = (T = eltype(x); T <: Real ? T : Bool)

# Distributions evaluate flat batches of array variates (the trailing
# dimensions are batch dimensions) directly where all variates lie in the
# support, univariate wrappers broadcast their point kernels:
for (bhead, phead) in ((:batched_logdensityof_impl, :logdensityof_impl), (:batched_logdensity_def, :logdensity_def))
    @eval function MeasureBase.$bhead(m::DistributionMeasure{<:ArrayLikeVariate{N}}, X::AbstractArray) where {N}
        _masked_batch_logpdf(m, X, Val(N))
    end
    @eval function MeasureBase.$bhead(m::DistributionMeasure{<:ArrayLikeVariate{0}}, X::AbstractArray)
        MeasureBase._scalar_kernel_broadcast(MeasureBase.$phead, m, X)
    end
end

function _masked_batch_logpdf(m::DistributionMeasure, X::AbstractArray, ::Val{N}) where {N}
    ndims(X) == N && return MeasureBase.logdensityof_impl(m, X)
    xs = eachslice(X; dims = ntuple(i -> N + i, ndims(X) - N))
    ins = map(x -> _logpdf_domain(m, x), xs)
    if all(ins)
        Distributions.logpdf(m.obj, X)
    else
        R = _logd_type(m.obj, X)
        map((x, i) -> i ? convert(R, Distributions.logpdf(m.obj, x)) : R(-Inf), xs, ins)
    end
end

@inline MeasureBase.unsafe_logdensityof(m::DistributionMeasure, x) = DensityInterface.logdensityof(m.obj, x)
@inline MeasureBase.insupportof(m::DistributionMeasure, x) = _dist_insupportof(m, m.obj, x)

# The support is unknown for distributions without an `insupport` method,
# decided at compile time (see `MeasureBase._has_batched_logdensity_def`):
@inline function _dist_insupportof(μ, d::D, x::X) where {D<:Distribution,X}
    if _has_insupport(D, X)
        Distributions.insupport(d, x) & _finite_variate(d, x)
    else
        NoFastInsupport{typeof(μ)}()
    end
end
@inline _has_insupport(::Type{D}, ::Type{X}) where {D,X} = @_adignore Core._hasmethod(Distributions.insupport, Tuple{D,X})

# Where densities are evaluated:
@inline _logpdf_domain(m::DistributionMeasure, x) = _known_or_finite(MeasureBase.insupportof(m, x), m.obj, x)
@inline _known_or_finite(ins, ::Distribution, x) = ins
@inline _known_or_finite(::NoFastInsupport, d::Distribution, x) = _finite_variate(d, x)

# Infinite values lie outside the support, where Distributions may evaluate
# to NaN:
@inline _finite_variate(::Distribution{Univariate}, x) = isfinite(x)
@inline _finite_variate(::Distribution{<:ArrayLikeVariate}, x) = all(isfinite, x)
@inline _finite_variate(::Distribution, x) = true

@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate{0},<:Continuous}) = Lebesgue()
@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate,<:Continuous}) = Lebesgue()^size(m.obj)
@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate{0},<:Discrete}) = Counting()
@inline MeasureBase.rootmeasure(m::DistributionMeasure{<:ArrayLikeVariate,<:Discrete}) = Counting()^size(m.obj)

# Distributions of non-array variates have no root measure (yet):
@inline MeasureBase.rootmeasure(m::DistributionMeasure) = _throw_no_rootmeasure(m.obj)
@noinline _throw_no_rootmeasure(d::Distribution) = throw(ArgumentError(
    "Wrapped distributions of $(nameof(Distributions.variate_form(typeof(d)))) variates have no root measure"
))

# Distributions on manifolds have densities relative to measures on their
# coordinates (e.g. Lebesgue measure on the free entries of symmetric
# matrices, or the surface measure of the sphere), which aren't declared as
# charts yet. Their opaque roots relate only to equal roots, so that they
# share the root with the distributions on the same manifold. Their degrees
# of freedom are the dimensions of the manifolds:
const _ManifoldDist = Union{
    Distributions.LKJ,Distributions.Wishart,Distributions.InverseWishart,
    Distributions.MatrixBeta,Distributions.MatrixFDist,Distributions.VonMisesFisher,
}

@inline function MeasureBase.rootmeasure(m::AsMeasure{<:_ManifoldDist})
    space, sz = _manifold(m.obj), size(m.obj)
    MeasureBase.OpaqueRoot(space, sz, _manifold_dim(space, sz))
end
@inline MeasureBase.getdof(m::AsMeasure{<:_ManifoldDist}) = getdof(rootmeasure(m))
@inline MeasureBase.fast_dof(m::AsMeasure{<:_ManifoldDist}) = MeasureBase.getdof(m)

# Correlation matrices, positive semi-definite matrices of a given rank, the
# unit sphere:
struct _CorrelationMatrices end
struct _PSDMatrices
    rank::Int
end
struct _UnitSphere end

@inline _manifold(::Distributions.LKJ) = _CorrelationMatrices()
@inline _manifold(d::Distributions.Wishart) = _PSDMatrices(LinearAlgebra.rank(d))
@inline _manifold(d::Union{Distributions.InverseWishart,Distributions.MatrixBeta,Distributions.MatrixFDist}) = _PSDMatrices(first(size(d)))
@inline _manifold(::Distributions.VonMisesFisher) = _UnitSphere()

@inline _manifold_dim(::_CorrelationMatrices, (n, _)::Dims{2}) = n * (n - 1) ÷ 2
@inline _manifold_dim(s::_PSDMatrices, (n, _)::Dims{2}) = n * s.rank - s.rank * (s.rank - 1) ÷ 2
@inline _manifold_dim(::_UnitSphere, (n,)::Dims{1}) = n - 1

# Logit-normal distributions on the simplex, like Dirichlet:
@inline MeasureBase.rootmeasure(m::AsMeasure{<:Distributions.MvLogitNormal}) = Lebesgue(MeasureBase.StdSimplex(length(m.obj)))
@inline MeasureBase.getdof(m::AsMeasure{<:Distributions.MvLogitNormal}) = length(m.obj) - 1
@inline MeasureBase.fast_dof(m::AsMeasure{<:Distributions.MvLogitNormal}) = length(m.obj) - 1

# Discrete distributions over integers count on all integers, those over
# other values on their atoms, so that the base measures of superpositions
# see these:
@inline MeasureBase.rootmeasure(m::AsMeasure{<:Distributions.DiscreteNonParametric}) = _atoms_root(m.obj)
@inline _atoms_root(::Distributions.DiscreteNonParametric{<:Integer}) = Counting()
@inline _atoms_root(d::Distributions.DiscreteNonParametric) = Counting(Distributions.support(d))

@inline MeasureBase.basemeasure(m::DistributionMeasure) = rootmeasure(m)

@inline MeasureBase.massof(::DistributionMeasure) = static(1.0)

@inline MeasureBase.mspace_elsize(d::Distribution) = MeasureBase.NoMSpaceElementSize{typeof(d)}()
@inline MeasureBase.mspace_flatsize(d::Distribution) = MeasureBase.NoMSpaceElementSize{typeof(d)}()
@inline MeasureBase.mspace_elsize(d::Distribution{Univariate}) = ()
@inline MeasureBase.mspace_elsize(d::Distribution{<:ArrayLikeVariate}) = size(d)
@inline MeasureBase.mspace_flatsize(d::Distribution{Univariate}) = ()
@inline MeasureBase.mspace_flatsize(d::Distribution{<:ArrayLikeVariate}) = size(d)
@inline MeasureBase.mspace_elsize(m::DistributionMeasure) = MeasureBase.mspace_elsize(m.obj)
@inline MeasureBase.mspace_flatsize(m::DistributionMeasure) = MeasureBase.mspace_flatsize(m.obj)
@inline MeasureBase.mspace_flatsize(::Type{<:Distribution{Univariate}}) = ()
@inline MeasureBase.mspace_ndims(::Type{<:Distribution{<:ArrayLikeVariate{N}}}) where {N} = N
@inline MeasureBase.mspace_ndims(::Type{AsMeasure{D}}) where {D<:Distribution} = MeasureBase.mspace_ndims(D)
@inline MeasureBase.mspace_flatsize(::Type{AsMeasure{D}}) where {D<:Distribution} = MeasureBase.mspace_flatsize(D)

@inline MeasureBase.preferred_stdmeasure(::Type{AsMeasure{D}}) where {D<:Distribution} = MeasureBase.preferred_stdmeasure(D)

@inline MeasureBase.getdof(m::DistributionMeasure{<:ArrayLikeVariate{0}}) = static(1)
@inline MeasureBase.getdof(m::DistributionMeasure) = _dist_getdof(m, Distributions.variate_form(typeof(m.obj)))
@inline _dist_getdof(m::DistributionMeasure, ::Type{<:ArrayLikeVariate}) = MeasureBase._default_getdof(typeof(m), basemeasure(m))
@inline _dist_getdof(m::DistributionMeasure, ::Type) = MeasureBase.NoDOF{typeof(m)}()
@inline MeasureBase.fast_dof(m::DistributionMeasure) = _dist_fast_dof(m, Distributions.variate_form(typeof(m.obj)))
@inline _dist_fast_dof(m::DistributionMeasure, ::Type{<:ArrayLikeVariate}) = MeasureBase._default_fast_dof(typeof(m), basemeasure(m))
@inline _dist_fast_dof(m::DistributionMeasure, ::Type) = MeasureBase.NoDOF{typeof(m)}()
@inline MeasureBase.fast_dof(m::DistributionMeasure{<:ArrayLikeVariate{0}}) = static(1)

# Delegate transport to the wrapped distribution:
@inline MeasureBase.transport_to_std(::Type{S}, m::DistributionMeasure, x) where {S<:StdMeasure} =
    MeasureBase.transport_to_std(S, m.obj, x)
@inline MeasureBase.transport_from_std(::Type{S}, m::DistributionMeasure, z) where {S<:StdMeasure} =
    MeasureBase.transport_from_std(S, m.obj, z)
@inline MeasureBase.batched_transport_to_std(::Type{S}, m::DistributionMeasure, X::AbstractArray) where {S<:StdMeasure} =
    MeasureBase.batched_transport_to_std(S, m.obj, X)
@inline MeasureBase.batched_transport_from_std(::Type{S}, m::DistributionMeasure, Z::AbstractArray) where {S<:StdMeasure} =
    MeasureBase.batched_transport_from_std(S, m.obj, Z)

@inline MeasureBase.paramnames(m::DistributionMeasure) = MeasureBase.paramnames(m.obj)
@inline MeasureBase.paramsof(m::DistributionMeasure) = MeasureBase.paramsof(m.obj)
