const EmptyNamedTuple = NamedTuple{(),Tuple{}}

function Base.show(io::IO, μ::AbstractMeasure)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)


export testvalue
testvalue(μ::AbstractMeasure) = testvalue(basemeasure(μ))
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

@generated function _rootmeasure(μ, ::StaticInt{n}) where {n}
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


@generated function instance(::Type{T}) where {T}
    return getfield(T, :instance)::T
end

# See https://github.com/cscherrer/KeywordCalls.jl/issues/22
@inline instance_type(f::F) where {F} = F
@inline instance_type(T::UnionAll) = Type{T}
@inline instance_type(T::DataType) = Type{T}

export basemeasure_depth

@inline function basemeasure_depth(μ::M) where {M}
    return basemeasure_depth(μ, basemeasure(μ))
end

@inline function basemeasure_depth(μ::M, β::M) where {M}
    return 0
end

@inline function basemeasure_depth(μ::M, β::B) where {M,B}
    return 1 + basemeasure_depth(β, basemeasure(β))
end

# Adapted from https://github.com/JuliaArrays/MappedArrays.jl/blob/46bf47f3388d011419fe43404214c1b7a44a49cc/src/MappedArrays.jl#L229
function func_string(f, types)
    ft = typeof(f)
    mt = ft.name.mt
    name = string(mt.name)
    if startswith(name, '#')
        # This is an anonymous function. See if it can be printed nicely
        lwrds = code_lowered(f, types)
        if length(lwrds) != 1
            return string(f)
        end
        lwrd = lwrds[1]
        c = lwrd.code
        if length(c) == 2 && ((isa(c[2], Expr) && c[2].head === :return) || (isdefined(Core, :ReturnNode) && isa(c[2], Core.ReturnNode)))
            # This is a single-line anonymous function, we should handle it
            s = lwrd.slotnames[2:end]
            result = ""
            if length(s) == 1
                result *= string(s[1])
                result *= "->"
            else
                result *= string(tuple(s...))
                result *= "->"
            end
            c1 = string(c[1])
            for i = 1:length(s)
                c1 = replace(c1, "_"*string(i+1)=>string(s[i]))
            end
            result *= c1
            return result
        else
            return string(f)
        end
    else
        return string(f)
    end
end

_eltype(T) = eltype(T)

function _eltype(g::Base.Generator{I}) where {I}
    Core.Compiler.return_type(g.f, Tuple{_eltype(I)})
end

function _eltype(::Type{Base.Generator{I,ComposedFunction{Outer,Inner}}}) where {Outer,Inner,I}
    _eltype(Base.Generator{_eltype(Base.Generator{I,Inner}), Outer})
end

function _eltype(::Type{Base.Generator{I,F}}) where {F<:Function,I}
    f = instance(F)
    Core.Compiler.return_type(f, Tuple{_eltype(I)})
end

function _eltype(::Type{Z}) where {Z<:Iterators.Zip}
    map(_eltype, Z.types[1].types)
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
