# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Samplers for discrete distributions on compute units other than the host
# and in compiled code. Host contexts draw with Distributions' samplers.
#
# The kernels are scalar functions of parameters and of uniform variates in
# `[0, 1)`, broadcast over batches, so parameters can be numbers or arrays
# on the compute unit. They compute in the float type of their arguments,
# the precision of the generative context, and return integer counts.
# Decimal constants are written as integer ratios, so that they don't
# promote their arguments.

# Means below `_small_mean` are sampled by inversion, larger ones by
# transformed rejection:
const _small_mean = 10

# The first `_inversion_terms` counts of a Poisson distribution with a rate
# of at most `_small_mean` have all but 1e-28 of the probability.
const _inversion_terms = 64

# Counts are `Int`s, and the acceptance tests take them in the float type
# `T` of the context, which holds integers exactly up to `maxintfloat(T)`:
@inline _max_count(::Type{T}) where {T<:AbstractFloat} = min(maxintfloat(T), T(typemax(Int)))

# Poisson rates stay 16 standard deviations below the largest count:
@inline _poisson_max_rate(::Type{T}) where {T<:AbstractFloat} = _max_count(T) - 16 * sqrt(_max_count(T))


# Inversion of the uniform variate `u` for a Poisson rate of at most
# `_small_mean`: the number of counts `k >= 1` with `cdf(k - 1) <= u`.
# Counts with a probability below `pfloor`, half the spacing of the values
# of `u`, are not reached: the largest `u`, in the last step of the cdf
# that the sum doesn't resolve, gets the last count with a resolvable
# probability.
@inline function _poisson_inversion(λ, u, pfloor)
    p = exp(-λ)
    s = p
    x = 0
    for k in 1:(_inversion_terms - 1)
        p = p * λ / k
        x += ifelse((u >= s) & (p >= pfloor), 1, 0)
        s += p
    end
    return x
end


# The transformed rejection method PTRS of W. Hörmann, "The transformed
# rejection method for generating Poisson random variables", Insurance:
# Mathematics and Economics 12 (1993), for rates of at least 10. A trial
# with the uniform variates `u` and `v` proposes the count
# `floor((2a / us + b) (u - 1/2) + λ + 0.43)` with `us = 1/2 - |u - 1/2|`
# and accepts it if `v invα / (a / us^2 + b)` is at most its probability.
#
# With the constants of the algorithm as it is commonly implemented (e.g.
# in NumPy), the envelope `invα / (a / us^2 + b)` lies up to 0.58 % below
# the probabilities of some counts for rates up to about 1e4 (at the upper
# ends of their intervals in `u`, 4e-7 of the probability mass at rate 10),
# so it is raised by 0.8 %. The squeeze and the tail test of the algorithm
# are left out: they only save evaluations of the probability, which
# batched code evaluates for every lane anyway.
@inline function _ptrs_constants(λ)
    b = (931 + 2530 * sqrt(λ)) / 1000
    a = (-5900 + 2483 * b) / 100000
    invα = 126 * (11239 + 113280 / (10 * b - 34)) / 1250000
    return b, a, invα
end

# The candidate count of a PTRS trial with the uniform variate `u`, -1 where
# the trial has none. The offset from the integer part of the rate is
# formed first and added in integers, so that candidates keep their
# resolution for large rates. Candidates above the largest count `kmax`
# count as none:
@inline function _ptrs_candidate(λ, u, kmax)
    b, a, _ = _ptrs_constants(λ)
    w = u - one(u) / 2
    us = one(u) / 2 - abs(w)
    us_safe = ifelse(us > 0, us, one(us))
    λ0 = floor(λ)
    d = floor((2 * a / us_safe + b) * w + (λ - λ0) + 43 * one(u) / 100)
    valid = (us > 0) & (d >= -λ0) & (d <= kmax - λ0)
    k = unsafe_trunc(Int, λ0) + unsafe_trunc(Int, ifelse(valid, d, zero(d)))
    return ifelse(valid, k, -1)
end

# Whether a PTRS trial with the uniform variates `u` and `v` accepts its
# candidate `k`:
@inline function _ptrs_accept(λ, u, v, k)
    b, a, invα = _ptrs_constants(λ)
    us = one(u) / 2 - abs(u - one(u) / 2)
    us_safe = ifelse(us > 0, us, one(us))
    v_safe = ifelse(v > 0, v, one(v))
    has_k = k >= 0
    kf = ifelse(has_k, k, 0) + zero(λ)
    return has_k & ((v <= 0) | (log(v_safe) + log(invα) - log(a / us_safe^2 + b) <= _poisson_logd(λ, kf)))
end

function _ptrs_round(ctx::GenContext, λ, dims::Dims)
    u = MeasureBase._rand_bulk(ctx, dims)
    v = MeasureBase._rand_bulk(ctx, dims)
    k = _ptrs_candidate.(λ, u, _max_count(get_precision(ctx)))
    return k, _ptrs_accept.(λ, u, v, k)
end


# Parameters take the precision of the context, host arrays of them go to
# its compute unit:
@inline _sampler_param(ctx::GenContext, x::Number) = convert_realtype(get_precision(ctx), x)
@inline _sampler_param(ctx::GenContext, x::AbstractArray) =
    MeasureBase._on_compute_unit(ctx, convert_realtype(get_precision(ctx), x))

@noinline _throw_poisson_rate(::Type{T}) where {T} = throw(ArgumentError(
    "The Poisson sampler for compute units other than the host supports rates from 0 to $(_poisson_max_rate(T)) in $T precision"
))

@inline _supported_rate(::Type{T}, λ) where {T} = (λ >= 0) & (λ <= _poisson_max_rate(T))

# Poisson variates of batch size `dims` for the rate or the array of rates
# `rate`, which broadcasts against the batch. Rates outside of the supported
# range throw, in compiled code their variates are -1.
function _poisson_draw(ctx::GenContext, rate, dims::Dims)
    T = get_precision(ctx)
    _poisson_draw_impl(ctx, T, _sampler_param(ctx, rate), dims)
end

function _poisson_draw_impl(ctx::GenContext, ::Type{T}, λ::AbstractFloat, dims::Dims) where {T}
    _supported_rate(T, λ) || _throw_poisson_rate(T)
    if λ < _small_mean
        _poisson_inversion.(λ, MeasureBase._rand_bulk(ctx, dims), eps(T) / 4)
    else
        MeasureBase._rejection_draw((c, p) -> _ptrs_round(c, p, dims), ctx, true, λ)
    end
end

function _poisson_draw_impl(ctx::GenContext, ::Type{T}, rate, dims::Dims) where {T}
    supported = _supported_rate.(T, rate)
    MeasureBase._values_known(rate) && !all(supported) && _throw_poisson_rate(T)
    λ = ifelse.(supported, rate, zero(T))
    small = λ .< _small_mean
    x_small = _poisson_inversion.(ifelse.(small, λ, zero(T)), MeasureBase._rand_bulk(ctx, dims), eps(T) / 4)
    λ_large = ifelse.(small, T(_small_mean), λ)
    x_large = MeasureBase._rejection_draw((c, p) -> _ptrs_round(c, p, dims), ctx, .!small, λ_large)
    return ifelse.(supported, ifelse.(small, x_small, x_large), -1)
end


function _nonnative_rand(ctx::GenContext, m::AsMeasure{<:Poisson}, sz::SizeLike)
    _poisson_draw(ctx, Distributions.rate(m.obj), MeasureBase.size_tuple(asnonstatic(sz)))
end

# Array products of Poisson distributions draw in one batch, with the rates
# of the marginals as an array of parameters:
const _PoissonProduct = MeasureBase.ProductMeasure{<:AbstractArray{<:AsMeasure{<:Poisson}}}

@inline _poisson_rates(mar::StructArray) = mar.obj.λ
@inline _poisson_rates(mar::AbstractArray) = map(m -> Distributions.rate(m.obj), mar)

function _poisson_product_draw(ctx::GenContext, d::_PoissonProduct, sz::SizeLike)
    rates = _poisson_rates(MeasureBase.marginals(d))
    batch = MeasureBase.size_tuple(asnonstatic(sz))
    _poisson_draw(ctx, reshape(rates, size(rates)..., map(_ -> 1, batch)...), (size(rates)..., batch...))
end

# Host contexts draw like other array products:
const _ArrayProduct = MeasureBase.ProductMeasure{<:AbstractArray}

function MeasureBase.rand_impl(ctx::GenContext, d::_PoissonProduct)
    MeasureBase._on_host(ctx) && return invoke(MeasureBase.rand_impl, Tuple{GenContext,_ArrayProduct}, ctx, d)
    _poisson_product_draw(ctx, d, ())
end

function MeasureBase.batched_rand_impl(ctx::GenContext, d::_PoissonProduct, sz::SizeLike)
    MeasureBase._on_host(ctx) &&
        return invoke(MeasureBase.batched_rand_impl, Tuple{GenContext,_ArrayProduct,SizeLike}, ctx, d, sz)
    _poisson_product_draw(ctx, d, sz)
end
