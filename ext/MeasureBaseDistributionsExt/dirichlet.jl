# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

const DirichletMeasure = AsMeasure{<:Dirichlet}

MeasureBase.getdof(d::Dirichlet) = length(d) - 1
MeasureBase.getdof(m::DirichletMeasure) = getdof(m.obj)
MeasureBase.fast_dof(m::DirichletMeasure) = getdof(m.obj)

@inline MeasureBase.preferred_stdmeasure(::Type{<:Dirichlet}) = StdLogistic

# Densities are relative to the coordinate Lebesgue measure on the simplex:
@inline MeasureBase.rootmeasure(m::DirichletMeasure) = Lebesgue(MeasureBase.StdSimplex(length(m.obj)))

# Test values are the barycentre (Distributions' samplers give NaN with
# constant draws):
@inline MeasureBase.testvalue(::Type{T}, m::DirichletMeasure) where {T} = testvalue(T, rootmeasure(m))

# Host compute units draw in log space, in the precision of the context or
# of the parameters, whichever is higher: log gamma variates, for
# concentrations `α < 1` with the boost `log(G′) - E / α` of a Gamma(1 + α)
# variate `G′` (`E` standard exponential), are normalized by a softmax. This
# is exact and avoids the underflow of the gamma variates themselves, which
# makes Distributions' sampler give NaN for small concentrations. Other
# compute units transport standard variates.
@inline MeasureBase.rand_impl(ctx::GenContext, m::DirichletMeasure) = MeasureBase.batched_rand_impl(ctx, m, ())

function MeasureBase.batched_rand_impl(ctx::GenContext, m::DirichletMeasure, sz::SizeLike)
    MeasureBase._on_host(ctx) || return MeasureBase._rand_default(ctx, m, sz, MeasureBase._NoRandImpl())
    T = float(promote_type(get_precision(ctx), Distributions.partype(m.obj)))
    α = m.obj.alpha
    X = Array{T}(undef, length(α), asnonstatic(sz)...)
    foreach(x -> _dirichlet_rand!(x, get_rng(ctx), α), eachslice(X; dims = ntuple(i -> i + 1, length(sz))))
    convert_realtype(get_precision(ctx), X)
end

function _dirichlet_rand!(x::AbstractVector{T}, rng::AbstractRNG, α::AbstractVector) where {T}
    # Where all `E / α` overflow, only the component with the smallest one
    # has a nonzero share:
    k, key_k = 0, T(Inf)
    for i in eachindex(x, α)
        a = convert(T, α[i])
        if a < 1
            e = randexp(rng, T)
            x[i] = _log_gamma_draw(rng, 1 + a) - e / a
            if x[i] == -Inf
                key = log(e) - log(a)
                key < key_k && ((k, key_k) = (i, key))
            end
        else
            x[i] = _log_gamma_draw(rng, a)
        end
    end
    m = maximum(x)
    if isfinite(m)
        x .= exp.(x .- m)
        x ./= sum(x)
    else
        x .= eachindex(x) .== k
    end
    return x
end

_dirichlet_rand!(x::AbstractVector{T}, ::ConstantRNG, ::AbstractVector) where {T} = fill!(x, inv(convert(T, length(x))))

# The log of a Gamma(a) variate for `a >= 1`, by Marsaglia & Tsang:
function _log_gamma_draw(rng::AbstractRNG, a::T) where {T<:AbstractFloat}
    d = a - one(T) / 3
    c = inv(sqrt(9 * d))
    while true
        x = randn(rng, T)
        v = (1 + c * x)^3
        v > 0 || continue
        u = rand(rng, T)
        (u < 1 - x^4 * 331 // 10000 || log(u) < x^2 / 2 + d - d * v + d * log(v)) && return log(d * v)
    end
end
