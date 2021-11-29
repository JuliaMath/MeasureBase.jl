import Base
using FillArrays: Fill
# """
# A power measure is a product of a measure with itself. The number of elements in
# the product determines the dimensionality of the resulting support.

# Note that power measures are only well-defined for integer powers.

# The nth power of a measure μ can be written μ^x.
# """
# PowerMeasure{M,N,D} = ProductMeasure{Fill{M,N,D}}

# function Base.show(io::IO, μ::PowerMeasure)
#     io = IOContext(io, :compact => true)
#     print(io, μ.data.value, " ^ ", size(μ.data))
# end

# function Base.show_unquoted(io::IO, μ::PowerMeasure{M,N,D}, indent::Int, prec::Int) where {M,N,D}
#     io = IOContext(io, :compact => true)
#     if Base.operator_precedence(:^) ≤ prec
#         print(io, "(")
#         show(io, μ.data.value)
#         print(io, ")")
#     else
#         show(io, size(μ.data))
#     end
#     return nothing
# end

export PowerMeasure

const PowerMeasure{M,N,R} = ProductMeasure{M,Kernel{Returns{M}, typeof(identity)}, LinearIndices{N, R}} 

Base.:^(μ::AbstractMeasure, ::Tuple{}) = μ

function Base.:^(μ::AbstractMeasure, dims::Integer...)
    return μ^dims
end

Base.show(io::IO, d::PowerMeasure) = print(io, d.f.f.value, " ^ ", size(d.xs))
Base.show(io::IO, d::PowerMeasure{M,1}) where {M} = print(io, d.f.f.value, " ^ ", length(d.xs))

# gentype(d::PowerMeasure{M,N}) where {M,N} = @inbounds Array{gentype(first(marginals(d))), N}

params(d::PowerMeasure) = params(first(marginals(d)))

# basemeasure(μ::PowerMeasure) = @inbounds basemeasure(first(μ.data))^size(μ.data)

@inline function basemeasure(d::PowerMeasure)
    _basemeasure(d, basemeasure(d.f.f.value))
end

@inline function _basemeasure(d::PowerMeasure, b)
    b ^ size(d.xs)
end

@inline function _basemeasure(d::PowerMeasure, b::FactoredBase)
    n = length(d.xs)
    inbounds(x) = all(b.inbounds, x)
    constℓ = n * b.constℓ
    varℓ() = n * b.varℓ()
    base = b.base^size(d.xs)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

function basemeasure_depth(d::PowerMeasure{M,N,R}) where {M,N,R}
    return basemeasure_depth(M)
end

function basemeasure_type(::Type{P}) where {M,N,R,P<:PowerMeasure{M,N,R}}
    return PowerMeasure{basemeasure_type(M), N, R}
end

function basemeasure_depth(::Type{P}) where {M<:PrimitiveMeasure,N,R,P<:PowerMeasure{M,N,R}}
    return static(0)
end

@inline function logdensity_def(d::PowerMeasure, x)
    d1 = d.f(first(d.xs))
    sum(xj -> logdensity_def(d1, xj), x)
end
