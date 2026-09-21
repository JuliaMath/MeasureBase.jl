@doc raw"""
    mkernel(f_β, f_c = OneTwoMany.secondarg)::Function

Constructs a generalized monadic transition kernel from a primary transition
kernel function `f_β` and a value combination function `f_c`.

`f_β` must behave like `β = f_β(a)`, taking a value `a` from a primary
measurable space and returning a measure-like object `β`.

`f_c` must behave like `c = f_c(a, b)`, taking a value `a` (like `f_β`) and
a value `b` from the measurable space of `β` and returning a value `c`.

`f_k = mkernel(f_β, f_c)` then acts like

```julia
f_k(a) ≡ pushfwd(c -> f_c(c[1], c[2]), productmeasure((Dirac(a), f_β(a))))
```

(`≡` denoting pseudocode-equivalency here). So with the default
`f_c == OneTwoMany.secondarg`, we just have `f_k(a) ≡ f_β(a)`.

Also,

```julia
mbind(mkernel(f_β, f_c), α) == mbind(f_β, α, f_c)
```

See also [`mbind`](@ref).
"""
function mkernel end
export mkernel


"""
    struct MeasureBase.MKernel <: Function

Represents a generalized monadic transition kernel.

User code should not create instances of `MKernel` directly, but should
call [`mkernel`](@ref) instead.
"""
struct MKernel{FT,FC} <: Function
    f_β::FT
    f_c::FC
end

(f_k::MKernel)(a) = mbind(f_k, Dirac(a))

@inline mkernel(f_β::MKernel) = f_β
@inline mkernel(f_β, f_c = secondarg) = _generic_mkernel_impl(f_β, f_c)

@inline _generic_mkernel_impl(f_β, f_c) = MKernel(f_β, f_c)
@inline _generic_mkernel_impl(f_β::MKernel, ::typeof(secondarg)) = f_β


@doc raw"""
    mbind(f_β, α::AbstractMeasure, f_c = OneTwoMany.secondarg)
    mbind(f_β::MeasureBase.MKernel, α::AbstractMeasure)

Constructs a monadic bind, resp. a hierarchical measure, from a transition
kernel function `f_β`, a primary measure `α` and a value combination
function `f_c`.

`f_β` must be a function that maps a point `a` from the space of the primary
measure `α` to a dependent secondary measure `β_a = f_β(a)`.
`ab = f_c(a, b)` must map such a point `a` and a point `b` from the
space of measure `β_a` to a combined value `ab = f_c(a, b)`.

The resulting measure

```julia
μ = mbind(f_β, α, f_c)
```

has the mathematical interpretation (on sets $$A$$ and $$B$$)

```math
\mu(f_c(A, B)) = \int_A \beta_a(B)\, \mathrm{d}\, \alpha(a)
```

When using the default `f_c = OneTwoMany.secondarg` (so `ab == b`) this
simplifies to

```math
\mu(B) = \int_A \beta_a(B)\, \mathrm{d}\, \alpha(a)
```

which is equivalent to a monadic bind, viewing measures as monads.

Computationally, `ab = rand(μ)` is equivalent to

```julia
a = rand(α)
β_a = f_β(a)
b = rand(β_a)
ab = f_c(a, b)
```

The measure `α` that went into the bind can be retrieved via
`boundmeasure(mbind(f_β, α, f_c)) == α` and the kernel via
`bindkernel(mbind(f_β, α, f_c)) == mkernel(f_β, f_c)`.

Densities on hierarchical measures can only be evaluated if `ab = f_c(a, b)`
can be unambiguously split into `a` and `b` again, knowing `α`. This is
currently implemented for `f_c` that is either `tuple` or `=>`/`Pair` (these
work for any combination of variate types), `vcat` (for tuple- or
vector-like variates) and `merge` (`NamedTuple` variates).
[`MeasureBase.tpmeasure_split_combined`](@ref) can be specialized to
support other choices for `f_c`.

# Extended help

Bayesian example with a correlated prior: Mathematically, let

    position = a1 ~ StdNormal()
    noise = a2 ~ pushforward(h(a1, ·), StdExponential())

where `h(a1, a2) = √(abs(a1) * a2)`. Because this prior on the space of
`A = A1 × A2 = (position, noise)` is a hierarchical measure (`a2` depends
on `a1`), we can construct it using `mbind` with `merge` as `f_c`:

```julia
using MeasureBase, AffineMaps

prior = mbind(
    productmeasure((
        position = StdNormal(),
    )), merge
) do a
    productmeasure((
        noise = pushfwd(setinverse(sqrt, setladj(x -> x^2, x -> log(2))) ∘ Mul(abs(a.position)), StdExponential()),
    ))
end

model = θ -> pushfwd(MulAdd(θ.noise, θ.position), StdNormal())^10

joint_θ_obs = mbind(model, prior, tuple)
prior_predictive = mbind(model, prior)

observation = rand(prior_predictive)
likelihood = likelihoodof(model, observation)

posterior = mintegrate(likelihood, prior)

θ = rand(prior)
logdensityof(posterior, θ)
```
"""
function mbind end
export mbind

@inline mbind(f_β) = Base.Fix1(mbind, f_β)

@inline function mbind(f_β, α::AbstractMeasure, f_c = secondarg)
    _generic_mbind_impl(f_β, asmeasure(α), f_c)
end

@inline function _generic_mbind_impl(f_β, α::AbstractMeasure, f_c)
    F, M, G = Core.Typeof(f_β), Core.Typeof(α), Core.Typeof(f_c)
    Bind{F,M,G}(f_β, α, f_c)
end

@inline _generic_mbind_impl(f_β, α::Dirac, f_c) = mcombine(f_c, α, asmeasure(f_β(α.x)))

@inline _generic_mbind_impl(@nospecialize(f_β), α::AbstractMeasure, ::typeof(firstarg)) = α
@inline _generic_mbind_impl(@nospecialize(f_β), α::Dirac, ::typeof(firstarg)) = α

@inline _generic_mbind_impl(f_k::MKernel, α::AbstractMeasure, ::typeof(secondarg)) =
    mbind(f_k.f_β, α, f_k.f_c)
@inline _generic_mbind_impl(f_k::MKernel, α::Dirac, ::typeof(secondarg)) =
    mbind(f_k.f_β, α, f_k.f_c)


"""
    struct MeasureBase.Bind <: AbstractMeasure

Represents a monadic bind resp. a hierarchical measure in general.

User code should not create instances of `Bind` directly, but should call
[`mbind`](@ref) instead.
"""
struct Bind{FT,M<:AbstractMeasure,FC} <: AbstractMeasure
    f_β::FT
    α::M
    f_c::FC
end

# ToDo: Store MKernel in Bind instead of separate fields f_β and f_c?


"""
    bindkernel(μ::Bind)::MKernel

Returns the monadic transition kernel of a monadic bind, so that
`bindkernel(mbind(f_k::MKernel, α)) == f_k`.

See [`mbind`](@ref) and [`mkernel`](@ref) for details.
"""
function bindkernel end
export bindkernel

bindkernel(μ::Bind) = mkernel(μ.f_β, μ.f_c)


"""
    boundmeasure(μ::Bind)::AbstractMeasure

Returns the measure that went into a monadic bind, so that
`boundmeasure(mbind(f_k, α)) == α`.

See [`mbind`](@ref) and [`mkernel`](@ref) for details.
"""
function boundmeasure end
export boundmeasure

boundmeasure(μ::Bind) = μ.α


_get_β_a(μ::Bind, a) = asmeasure(μ.f_β(a))

function transportmeasure(μ::Bind, x)
    tpm_α, a, b = tpmeasure_split_combined(μ.f_c, μ.α, x)
    tpm_β_a = transportmeasure(_get_β_a(μ, a), b)
    mcombine(μ.f_c, tpm_α, tpm_β_a)
end

localmeasure(μ::Bind, x) = transportmeasure(μ, x)

tpmeasure_split_combined(f_c, μ::Bind, xy) = _bind_tpm_sc(f_c, μ, xy)

function _bind_tpm_sc(::typeof(tuple), μ::Bind, xy::Tuple{Vararg{Any,2}})
    x, y = xy[1], xy[2]
    tpm_μ = transportmeasure(μ, x)
    return tpm_μ, x, y
end

function _bind_tpm_sc(::Type{Pair}, μ::Bind, xy::Pair)
    x, y = xy.first, xy.second
    tpm_μ = transportmeasure(μ, x)
    return tpm_μ, x, y
end

const _BindBy{FC} = Bind{<:Any,<:AbstractMeasure,FC}

@inline preferred_stdmeasure(::Type{<:Bind{<:Any,M}}) where {M} = preferred_stdmeasure(M)
_bind_tpm_sc(f_c::typeof(vcat), μ::_BindBy{typeof(vcat)}, xy::AbstractVector) =
    _bind_tpm_sc_cat(f_c, μ, xy)
_bind_tpm_sc(f_c::typeof(merge), μ::_BindBy{typeof(merge)}, xy::NamedTuple) =
    _bind_tpm_sc_cat(f_c, μ, xy)

function _bind_tpm_sc_cat_lμabyxy(f_c, μ, xy)
    tpm_α, a, by = tpmeasure_split_combined(μ.f_c, μ.α, xy)
    β_a = _get_β_a(μ, a)
    tpm_β_a, b, y = tpmeasure_split_combined(f_c, β_a, by)
    tpm_μ = mcombine(μ.f_c, tpm_α, tpm_β_a)
    return tpm_μ, a, b, y, xy
end

function _bind_tpm_sc_cat(f_c::typeof(vcat), μ::_BindBy{typeof(vcat)}, xy::AbstractVector)
    tpm_μ, a, b, y, xy = _bind_tpm_sc_cat_lμabyxy(f_c, μ, xy)
    # Don't use `x = f_c(a, b)` here, would allocate, splitting xy can use views:
    x, y = split_at(xy, maybestatic_length(a) + maybestatic_length(b))
    return tpm_μ, x, y
end

function _bind_tpm_sc_cat(f_c::typeof(merge), μ::_BindBy{typeof(merge)}, xy::NamedTuple)
    tpm_μ, a, b, y, xy = _bind_tpm_sc_cat_lμabyxy(f_c, μ, xy)
    return tpm_μ, f_c(a, b), y
end


@inline insupport(μ::Bind, ::Any) = NoFastInsupport{typeof(μ)}()

@inline getdof(μ::Bind) = NoDOF{typeof(μ)}()

# Bypass `checked_arg`, would require potentially costly evaluation of f_β:
@inline checked_arg(::Bind, x) = x

rootmeasure(::Bind) =
    throw(ArgumentError("root measure is implicit, but can't be instantiated, for Bind"))

basemeasure(::Bind) = throw(ArgumentError("basemeasure is not available for Bind"))

# Test values follow the primary test value through the secondary measure:
function testvalue(::Type{T}, μ::Bind) where {T}
    a = testvalue(T, μ.α)
    _combine_variates(μ.f_c, a, testvalue(T, _get_β_a(μ, a)))
end

logdensity_def(::Bind, x) =
    throw(ArgumentError("logdensity_def is not available for Bind"))

# Density evaluation consumes the variate parts of the primary and secondary
# measure in a single pass, using the with-rest protocol for value-dependent
# variate sizes:

logdensityof_impl(μ::Bind, x) = _bind_ld_impl(μ.f_c, μ, x)

unsafe_logdensityof(μ::Bind, x) = logdensityof_impl(μ, x)

function _bind_ld_impl(::typeof(tuple), μ::Bind, xy::Tuple{Vararg{Any,2}})
    a, b = xy[1], xy[2]
    logdensityof(μ.α, a) + logdensityof(_get_β_a(μ, a), b)
end

function _bind_ld_impl(::Type{Pair}, μ::Bind, xy::Pair)
    a, b = xy.first, xy.second
    logdensityof(μ.α, a) + logdensityof(_get_β_a(μ, a), b)
end

function _bind_ld_impl(::Union{typeof(vcat),typeof(merge)}, μ::Bind, xy)
    ℓ, _, x_rest = logdensityof_with_rest(μ, xy)
    isempty(x_rest) || _throw_stream_too_long()
    return ℓ
end

function _bind_ld_impl(@nospecialize(f_c), @nospecialize(μ::Bind), @nospecialize(xy))
    throw(
        ArgumentError(
            "Can't compute density of a bind with value combination function of type $(nameof(typeof(f_c)))",
        ),
    )
end

# The secondary measure depends on the primary variate, so streams are
# consumed one by one:
@inline fixed_stream_size(::Type{<:Bind}) = static(false)
@inline mspace_ndims(::Type{<:_BindBy{typeof(vcat)}}) = 1

function logdensityof_with_rest(μ::_BindBy{typeof(vcat)}, x::AbstractVector)
    ℓ_a, a, x2 = logdensityof_with_rest(μ.α, x)
    ℓ_b, b, x_rest = logdensityof_with_rest(_get_β_a(μ, a), x2)
    x_μ, _ = split_at(x, maybestatic_length(x) - maybestatic_length(x_rest))
    return ℓ_a + ℓ_b, x_μ, x_rest
end

function logdensityof_with_rest(μ::_BindBy{typeof(merge)}, x::NamedTuple)
    ℓ_a, a, x2 = logdensityof_with_rest(μ.α, x)
    ℓ_b, b, x_rest = logdensityof_with_rest(_get_β_a(μ, a), x2)
    return ℓ_a + ℓ_b, merge(a, b), x_rest
end

function batched_logdensityof_with_rest(μ::Bind, x::AbstractVector, ::Tuple{})
    ℓ, _, x_rest = logdensityof_with_rest(μ, x)
    return ℓ, x_rest
end

batched_logdensityof_impl(μ::_BindBy{typeof(vcat)}, X::AbstractArray) = _streamwise_ld(logdensityof_impl, μ, X)

# Batches of streams containing binds are consumed stream by stream (by
# the outermost stream combinator, see `fixed_stream_size`):
@noinline function batched_logdensityof_with_rest(::Bind, ::AbstractArray, ::SizeLike)
    throw(ArgumentError("Batches of variate streams containing binds must be consumed stream by stream"))
end
batched_logdensityof_impl(μ::_BindBy{typeof(vcat)}, x::AbstractVector) = _bind_ld_impl(vcat, μ, x)


function rand_impl(ctx::GenContext, μ::Bind)
    a = rand_impl(ctx, μ.α)
    b = rand_impl(ctx, _get_β_a(μ, a))
    return _combine_variates(μ.f_c, a, b)
end

# The secondary measure depends on the primary variate, so batches are
# generated variate by variate:
batched_rand_impl(ctx::GenContext, μ::Bind, sz::SizeLike) = _batched_rand_pointwise(ctx, μ, sz)


# Transport consumes the variate parts of the primary and secondary
# measure in a single pass, analogous to density evaluation:

transport_to_std(::Type{S}, μ::Bind, ab) where {S<:StdMeasure} = _bind_to_std(S, μ.f_c, μ, ab)

function _bind_to_std(::Type{S}, f_c, μ::Bind, ab) where {S}
    tpm_α, a, b = tpmeasure_split_combined(f_c, μ.α, ab)
    vcat(_as_stdstream(transport_to_std(S, tpm_α, a)), _as_stdstream(transport_to_std(S, _get_β_a(μ, a), b)))
end

function _bind_to_std(::Type{S}, ::Union{typeof(vcat),typeof(merge)}, μ::Bind, ab) where {S}
    z, _, x_rest = transport_to_std_with_rest(S, μ, ab)
    if !isempty(x_rest)
        throw(ArgumentError("Variate too long during transport of a bind"))
    end
    return z
end

function transport_to_std_with_rest(::Type{S}, μ::_BindBy{typeof(vcat)}, x::AbstractVector) where {S<:StdMeasure}
    z_a, a, x2 = transport_to_std_with_rest(S, μ.α, x)
    z_b, _, x_rest = transport_to_std_with_rest(S, _get_β_a(μ, a), x2)
    x_μ, _ = split_at(x, maybestatic_length(x) - maybestatic_length(x_rest))
    return vcat(z_a, z_b), x_μ, x_rest
end

function transport_to_std_with_rest(::Type{S}, μ::_BindBy{typeof(merge)}, x::NamedTuple) where {S<:StdMeasure}
    z_a, a, x2 = transport_to_std_with_rest(S, μ.α, x)
    z_b, b, x_rest = transport_to_std_with_rest(S, _get_β_a(μ, a), x2)
    return vcat(z_a, z_b), merge(a, b), x_rest
end

function transport_from_std_with_rest(::Type{S}, μ::Bind, z::AbstractVector) where {S<:StdMeasure}
    a, z2 = transport_from_std_with_rest(S, μ.α, z)
    b, z_rest = transport_from_std_with_rest(S, _get_β_a(μ, a), z2)
    return _combine_variates(μ.f_c, a, b), z_rest
end

function transport_from_std(::Type{S}, μ::Bind, z::AbstractVector) where {S<:StdMeasure}
    x, z_rest = transport_from_std_with_rest(S, μ, z)
    isempty(z_rest) || _throw_std_length_mismatch()
    return x
end

# The secondary measure depends on the primary variate, so batches of
# streams are consumed stream by stream (by the outermost stream
# combinator, see `fixed_stream_size`):
function batched_transport_to_std_with_rest(::Type{S}, μ::Bind, X::AbstractArray, sz::SizeLike) where {S<:StdMeasure}
    _bind_to_std_with_rest(S, μ, X, sz)
end
function _bind_to_std_with_rest(::Type{S}, μ::Bind, x::AbstractVector, ::Tuple{}) where {S}
    z, _, x_rest = transport_to_std_with_rest(S, μ, x)
    return z, x_rest
end
@noinline function _bind_to_std_with_rest(::Type{S}, ::Bind, ::AbstractArray, ::SizeLike) where {S}
    throw(ArgumentError("Batches of variate streams containing binds must be consumed stream by stream"))
end
