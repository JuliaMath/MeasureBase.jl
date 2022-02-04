
export For
using Random
import Base

struct For{T, F, I} <: AbstractProductMeasure
    f::F
    inds::I

    @inline function For{T}(f::F, inds::I) where {T,F,I<:Tuple}
        new{T,instance_type(f),I}(f, inds)
    end

    @inline For{T,F,I}(f::F, inds::I) where {T,F,I} = new{T,F,I}(f,inds)

    function For{Union{},F,I}(f::F, inds::I) where {F,I}
        println("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
        println.(stacktrace())
        println("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
        @show f(first(zip(inds...))...)
        @error "Empty `For` construction"
    end
end


@generated function For(f::F, inds::I) where {F,I<:Tuple}
    eltypes = Tuple{eltype.(I.types)...}
    quote
        $(Expr(:meta, :inline))
        T = Core.Compiler.return_type(f, $eltypes)
        For{T,F,I}(f, inds)
    end
end


# For(f, gen::Base.Generator) = ProductMeasure(Base.Generator(f ∘ gen.f, gen.iter))

@inline function logdensity_def(d::For{T,F,I}, x::AbstractVector{X}) where {X,T,F,I<:Tuple{<:AbstractVector}}
    ℓ = 0.0
    @inbounds for j in eachindex(x)
        ℓ += logdensity_def(d.f(j), x[j])
    end
    ℓ
end

function logdensity_def(d::For, x::AbstractVector) 
    sum(eachindex(x)) do i
        @inbounds logdensity_def(d.f(getindex.(d.inds,i)...), x[i])
    end
end

function logdensity_def(d::For{T,F,I}, x::AbstractArray{X}) where {T,F,I,X}
    ℓ = zero(typeintersect(AbstractFloat, Core.Compiler.return_type(logdensity_def, Tuple{T,X})))

    @inbounds for j in CartesianIndices(x)
        i = (getindex(ind, j) for ind in d.inds)
        ℓ += logdensity_def(d.f(i...), x[j])
    end
    ℓ
end

function logdensity_def(d::For{T,F,I}, x) where {N,T,F,I<:NTuple{N,<:Base.Generator}}
    sum(zip(x, d.inds...)) do (xⱼ, dⱼ...)
        logdensity_def(d.f(dⱼ...), xⱼ)
    end
end

function logdensity_def(d::For{T,F,I}, x::AbstractVector) where {N,T,F,I<:NTuple{N,<:Base.Generator}}
    sum(zip(x, d.inds...)) do (xⱼ, dⱼ...)
        logdensity_def(d.f(dⱼ...), xⱼ)
    end
end

function marginals(d::For{T,F,Tuple{I}}) where {T,F,I}
    f = d.f
    data = first(d.inds)
    MappedArrays.ReadonlyMappedArray{T,ndims(data),typeof(data),typeof(f)}(f, data)
end

function marginals(d::For{T,F,I}) where {T,F,I}
    f = d.f
    data = d.inds
    MappedArrays.ReadonlyMultiMappedArray{T,ndims(first(data)),typeof(data),typeof(f)}(f, data)
end

function marginals(d::For{T,F,I}) where {N,T,F,I<:NTuple{N,<:Base.Generator}}
    Iterators.map(d.f, d.inds...)
end

@inline function basemeasure(d::For{T,F,I}) where {T,F,I}
    B = Core.Compiler.return_type(basemeasure, Tuple{T})
    _basemeasure(d, B, static(Base.issingletontype(B)))
end

@inline function _basemeasure(d::For{T,F,I}, ::Type{B}, ::True) where {T,F,I,B}
    instance(B) ^ axes(d.inds)
end

@inline function _basemeasure(d::For{T,F,I}, ::Type{B}, ::False) where {T,F,I,B<:AbstractMeasure}
    f = basekleisli(d.f)
    For{B}(f, d.inds)
end

@inline function _basemeasure(d::For{T,F,I}, ::Type{B}, ::False) where {T,F,I,B}
    productmeasure(basemeasure.(marginals(d)))
end

function _basemeasure(d::For{T,F,I}, ::Type{B}, ::True) where {N,T<:AbstractMeasure,F,I<:NTuple{N,<:Base.Generator},B}
    return instance(B) ^ minimum(length, d.inds)
end

function _basemeasure(d::For{T,F,I}, ::Type{B}, ::False) where {N,T<:AbstractMeasure,F,I<:NTuple{N,<:Base.Generator},B}
    f = basekleisli(d.f)
    For{B}(f, d.inds)
end

function Pretty.tile(d::For{T}) where {T}
    result = Pretty.literal("For{")
    result *= Pretty.tile(T)
    result *= Pretty.literal("}")
    result *= Pretty.list_layout(
        [
            Pretty.literal(func_string(d.f, Tuple{_eltype.(d.inds)...})),
            Pretty.tile.(d.inds)...
        ]
    )
end



"""
    For(f, base...)

`For` provides a convenient way to construct a `ProductMeasure`. There are
several options for the `base`. With Julia's `do` notation, this can look very
similar to a standard `for` loop, while maintaining semantics structure that's
easier to work with.

------------

# `For(f, base::Int...)`

When one or several `Int` values are passed for `base`, the result is treated as
depending on `CartesianIndices(base)`. 

```
julia> For(3) do λ Exponential(λ) end |> marginals
3-element mappedarray(MeasureBase.var"#17#18"{var"#15#16"}(var"#15#16"()), ::CartesianIndices{1, Tuple{Base.OneTo{Int64}}}) with eltype Exponential{(:λ,), Tuple{Int64}}:
 Exponential(λ = 1,)
 Exponential(λ = 2,)
 Exponential(λ = 3,)
```

```
julia> For(4,3) do μ,σ Normal(μ,σ) end |> marginals
4×3 mappedarray(MeasureBase.var"#17#18"{var"#11#12"}(var"#11#12"()), ::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}}) with eltype Normal{(:μ, :σ), Tuple{Int64, Int64}}:
 Normal(μ = 1, σ = 1)  Normal(μ = 1, σ = 2)  Normal(μ = 1, σ = 3)
 Normal(μ = 2, σ = 1)  Normal(μ = 2, σ = 2)  Normal(μ = 2, σ = 3)
 Normal(μ = 3, σ = 1)  Normal(μ = 3, σ = 2)  Normal(μ = 3, σ = 3)
 Normal(μ = 4, σ = 1)  Normal(μ = 4, σ = 2)  Normal(μ = 4, σ = 3)
```

-------

# `For(f, base::AbstractArray...)``

In this case, `base` behaves as if the arrays are `zip`ped together before
applying the map.

```
julia> For(randn(3)) do x Exponential(x) end |> marginals
3-element mappedarray(x->Main.Exponential(x), ::Vector{Float64}) with eltype Exponential{(:λ,), Tuple{Float64}}:
 Exponential(λ = -0.268256,)
 Exponential(λ = 1.53044,)
 Exponential(λ = -1.08839,)
```

```
julia> For(1:3, 1:3) do μ,σ Normal(μ,σ) end |> marginals
3-element mappedarray((:μ, :σ)->Main.Normal(μ, σ), ::UnitRange{Int64}, ::UnitRange{Int64}) with eltype Normal{(:μ, :σ), Tuple{Int64, Int64}}:
 Normal(μ = 1, σ = 1)
 Normal(μ = 2, σ = 2)
 Normal(μ = 3, σ = 3)
```

----

# `For(f, base::Base.Generator)`

For `Generator`s, the function maps over the values of the generator:

```
julia> For(eachrow(rand(4,2))) do x Normal(x[1], x[2]) end |> marginals |> collect
4-element Vector{Normal{(:μ, :σ), Tuple{Float64, Float64}}}:
 Normal(μ = 0.255024, σ = 0.570142)
 Normal(μ = 0.970706, σ = 0.0776745)
 Normal(μ = 0.731491, σ = 0.505837)
 Normal(μ = 0.563112, σ = 0.98307)
```

"""
@inline For{T}(f, inds...) where {T} = For{T}(f, inds)
@inline For{T}(f, n::Integer) where {T} = For{T}(f, Base.OneTo(n))
@inline For{T}(f, inds::Integer...) where {T} = For{T}(i -> f(Tuple(i)...), Base.CartesianIndices(inds))

@inline For(f, inds...) = For(f, inds)
@inline For(f, n::Integer) = For(f, Base.OneTo(n))
@inline For(f, inds::Integer...) = For(i -> f(Tuple(i)...), Base.CartesianIndices(inds))
# For(f, inds::Base.Generator) = productmeasure(mymap(f, inds))

function Random.rand!(rng::AbstractRNG, d::For{T,F,I}, x) where {T,F,I} 
    mar = marginals(d)
    @inbounds for (dⱼ, j) in zip(mar, eachindex(x))
        x[j] = rand(rng,dⱼ)
    end
    return x
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::For{M,F,I}) where {T,M,F,I}
    _rand_product(rng, T, marginals(d), M)
end