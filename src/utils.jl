function Base.show(io::IO, μ::AbstractMeasure)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

export testvalue

@inline testvalue(μ) = testvalue(Float64, μ)

@inline testvalue(::Type{T}, μ) where {T} = rand(ConstantRNG(), T, μ)

testvalue(::Type{T}) where {T} = zero(T)

export rootmeasure

"""
    rootmeasure(μ::AbstractMeasure)

It's sometimes important to be able to find the fix point of a measure under
`basemeasure`. That is, to start with some measure and apply `basemeasure`
repeatedly until there's no change. That's what this does.
"""
@inline function rootmeasure(μ)
    n = basemeasure_depth(μ)
    _rootmeasure(μ, static(n))
end

@generated function _rootmeasure(μ, ::StaticInteger{n}) where {n}
    q = quote end
    foreach(1:n) do _
        push!(q.args, :(μ = basemeasure(μ)))
    end
    push!(q.args, :(return μ))
    return q
end

function instance(::Type{T}) where {T}
    return getfield(T, :instance)::T
end

export basemeasure_depth

@inline function basemeasure_depth(μ::M) where {M}
    b_0 = μ
    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        b_{i} = basemeasure(b_{i - 1})
        if b_{i} isa typeof(b_{i - 1})
            return static(i - 1)
        end
    end
    return static(10)
end

"""
    basemeasure_sequence(m)

Construct the longest `Tuple` starting with `m` having each term as the base
measure of the previous term, and with no repeated entries.
"""
@inline function basemeasure_sequence(μ::M) where {M}
    b_1 = μ
    done = false
    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        b_{i + 1} = if done
            nothing
        else
            basemeasure(b_{i})
        end
        if b_{i + 1} isa typeof(b_{i})
            done = true
            b_{i + 1} = nothing
        end
    end
    return filter(!isnothing, Base.Cartesian.@ntuple 10 b)
end

rmap(f, x) = f(x)

function rmap(f, t::Tuple)
    map(x -> rmap(f, x), t)
end

function rmap(f, nt::NamedTuple{N,T}) where {N,T}
    NamedTuple{N}(map(x -> rmap(f, x), values(nt)))
end

insupport(m::AbstractMeasure) = Base.Fix1(insupport, m)

using InverseFunctions: FunctionWithInverse

unwrap(f) = f
unwrap(f::FunctionWithInverse) = f.f

fcomp(f, g) = fchain(g, f)
fcomp(::typeof(identity), g) = g
fcomp(f, ::typeof(identity)) = f
fcomp(::typeof(identity), ::typeof(identity)) = identity

isneginf(x) = isinf(x) & (x < zero(x))
isposinf(x) = isinf(x) & (x > zero(x))

# `cond ? f_true() : f_false()` for plain booleans, so the value not selected
# is never computed. Other conditions (traced or on devices) select
# branch-free between both values:
@inline _select(cond::Bool, f_true::F, f_false::G) where {F,G} = cond ? f_true() : f_false()
@inline function _select(cond, f_true::F, f_false::G) where {F,G}
    y = f_true()
    ifelse(cond, y, oftype(y, f_false()))
end

containsnan(x::Number) = isnan(x)
containsnan(x) = any(containsnan, x)


# ForwardDiffPullbacks dummy `fwddiff`, overloaded by
# ForwardDiffPullbacks extension when loaded:
@inline _fwddiff(f) = f


# Autodiff ignore:

@inline _adignore_call(f) = f()

macro _adignore(expr)
    :(_adignore_call(() -> $(esc(expr))))
end


"""
    MeasureBase.convert_realtype(::Type{T}, x) where {T<:Real}

Convert the floating-point numbers in `x` to `T`. Other numbers, like
integers and dual numbers, are kept as they are.
"""
function convert_realtype end

@inline convert_realtype(::Type{T}, x::T) where {T<:Real} = x
@inline convert_realtype(::Type{T}, x::AbstractFloat) where {T<:Real} = T(x)
@inline convert_realtype(::Type{T}, x::Number) where {T<:Real} = x
@inline convert_realtype(::Type{T}, x::AbstractArray{T}) where {T<:Real} = x
convert_realtype(::Type{T}, x::AbstractArray{<:AbstractFloat}) where {T<:Real} = T.(x)
@inline convert_realtype(::Type{T}, x::AbstractArray{<:Number}) where {T<:Real} = x
convert_realtype(::Type{T}, x::Union{Tuple,NamedTuple}) where {T<:Real} =
    map(Base.Fix1(convert_realtype, T), x)
convert_realtype(::Type{T}, x::AbstractArray) where {T<:Real} =
    map(Base.Fix1(convert_realtype, T), x)

# Regularized incomplete gamma and beta functions and their inverses, with
# the log-densities of the standard gamma and beta distributions for their
# derivatives. Implemented in the Distributions extension, differentiated
# with respect to the variate resp. probability argument in the autodiff
# extensions:
function _gamma_cdf end
function _gamma_quantile end
function _beta_cdf end
function _beta_quantile end
function _gamma_logpdf end
function _beta_logpdf end
function _gamma_cdf_impl end
function _gamma_quantile_impl end
function _beta_cdf_impl end
function _beta_quantile_impl end

# The dual number type among the arguments of such a function, `Nothing`
# for plain numbers (the ForwardDiff extension adds dual numbers):
@inline _dualtag() = Nothing
@inline _dualtag(::Number, rest::Number...) = _dualtag(rest...)

# Distributions implementation hooks, specialized for dual numbers in the
# ForwardDiff extension:
function _trafo_logcdf_impl end
function _trafo_logccdf_impl end
function _trafo_quantile_impl end
function _trafo_cquantile_impl end
function _dist_quantile end
function _dist_cquantile end
function _dist_params_numtype end
