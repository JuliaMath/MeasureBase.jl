const EmptyNamedTuple = NamedTuple{(),Tuple{}}
const NonEmptyTuple = Tuple{Any,Vararg{Any}}

function Base.show(io::IO, μ::AbstractMeasure)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)

export testvalue

@inline testvalue(μ) = testvalue(Float64, μ)

@inline testvalue(::Type{T}, μ) where {T} = rand(ConstantRNG(), T, μ)

testvalue(::Type{T}) where {T} = zero(T)

export rootmeasure

basemeasure(μ, x) = basemeasure(μ)

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
    return q
end

# Base on the Tricks.jl README
using Tricks
struct Iterable end
struct NonIterable end
function isiterable(::Type{T}) where {T}
    static_hasmethod(iterate, Tuple{T}) ? Iterable() : NonIterable()
end

# issingletontype(@nospecialize(t)) = (@_pure_meta; isa(t, DataType) && isdefined(t, :instance))

# @generated function instance(::Type{T}) where {T}
#     return getfield(T, :instance)::T
# end

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

commonbase(μ, ν) = commonbase(μ, ν, Any)

"""
    commonbase(μ, ν, T) -> Tuple{StaticInt{i}, StaticInt{j}}

Find minimal (with respect to their sum) `i` and `j` such that there is a method

    logdensity_def(basemeasure_sequence(μ)[i], basemeasure_sequence(ν)[j], ::T)

This is used in `logdensity_rel` to help make that function efficient.
"""
@inline function commonbase(μ, ν, ::Type{T}) where {T}
    return commonbase(basemeasure_sequence(μ), basemeasure_sequence(ν), T)
end

@generated function commonbase(μ::M, ν::N, ::Type{T}) where {M<:Tuple,N<:Tuple,T}
    m = schema(M)
    n = schema(N)

    sols = Iterators.filter(
        ((i, j),) -> static_hasmethod(logdensity_def, Tuple{m[i],n[j],T}),
        Iterators.product(1:length(m), 1:length(n)),
    )
    isempty(sols) && return :(nothing)
    minsol = static.(argmin(((i, j),) -> i + j, sols))
    quote
        $minsol
    end
end

mymap(f, gen::Base.Generator) = mymap(f ∘ gen.f, gen.iter)
mymap(f, inds...) = Iterators.map(f, inds...)

function infer_zero(f, args...)
    inferred_type = Core.Compiler.return_type(f, typeof.(args))
    zero(typeintersect(AbstractFloat, inferred_type))
end

@inline function allequal(f, x::AbstractArray)
    val = f(first(x))
    @simd for xj in x
        f(xj) == val || return false
    end
    return true
end

allequal(x::AbstractArray) = allequal(identity, x)

rmap(f, x) = f(x)

function rmap(f, t::Tuple)
    map(x -> rmap(f, x), t)
end

function rmap(f, nt::NamedTuple{N,T}) where {N,T}
    NamedTuple{N}(map(x -> rmap(f, x), values(nt)))
end

insupport(m::AbstractMeasure) = Base.Fix1(insupport, m)

@inline return_type(f, args::Tuple) = Core.Compiler.return_type(f, Tuple{typeof.(args)...})

unstatic(::Type{T}) where {T} = T
unstatic(::Type{StaticFloat64{X}}) where {X} = Float64

using InverseFunctions: FunctionWithInverse

unwrap(f) = f
unwrap(f::FunctionWithInverse) = f.f

fcomp(f, g) = fchain(g, f)
fcomp(::typeof(identity), g) = g
fcomp(f, ::typeof(identity)) = f
fcomp(::typeof(identity), ::typeof(identity)) = identity

near_neg_inf(::Type{T}) where {T<:Real} = T(-1E38) # Still fits into Float32

isneginf(x) = isinf(x) && x < zero(x)
isposinf(x) = isinf(x) && x > zero(x)

isapproxzero(x::T) where {T<:Real} = x ≈ zero(T)
isapproxzero(A::AbstractArray) = all(isapproxzero, A)

isapproxone(x::T) where {T<:Real} = x ≈ one(T)
isapproxone(A::AbstractArray) = all(isapproxone, A)
