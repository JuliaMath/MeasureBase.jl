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

const PowerMeasure{F,S,T,N,A} = ProductMeasure{F,S,Fill{T,N,A}}

Base.:^(μ::AbstractMeasure, ::Tuple{}) = μ

function Base.:^(μ::AbstractMeasure, dims::Integer...)
    return μ ^ dims
end

function Base.:^(μ::M, dims::NTuple{N,I}) where {M<:AbstractMeasure,N,I<:Integer}
    powermeasure(μ, dims)
end

# Same as PowerMeasure
function Base.show(io::IO, d::ProductMeasure{F,S,<:Fill}) where {F,S}
    io = IOContext(io, :compact => true)
    print(io, d.f(first(d.pars)), " ^ ", size(d.pars))
end

# Same as PowerMeasure{F,T,1} where {F,T}
function Base.show(io::IO, d::ProductMeasure{F,S,Fill{T,1,A}}) where {F,S,T,A}
    io = IOContext(io, :compact => true)
    print(io, d.f(first(d.pars)), " ^ ", size(d.pars)[1])
end

# sampletype(d::PowerMeasure{M,N}) where {M,N} = @inbounds Array{sampletype(first(marginals(d))), N}



params(d::ProductMeasure{F,S,<:Fill}) where {F,S} = params(first(marginals(d)))

params(::Type{P}) where {F,S,P<:ProductMeasure{F,S,<:Fill}} = params(D)

# basemeasure(μ::PowerMeasure) = @inbounds basemeasure(first(μ.data))^size(μ.data)

# Same as PowerMeasure
@inline function basemeasure(d::ProductMeasure{F,S,<:Fill}) where {F,S}
    _basemeasure(d, (basemeasure(d.f(first(d.pars)))))
end

# Same as PowerMeasure
@inline function _basemeasure(d::ProductMeasure{F,S,<:Fill}, b) where {F,S}
    b ^ size(d.pars)
end

# Same as PowerMeasure
@inline function _basemeasure(d::ProductMeasure{F,S,<:Fill}, b::FactoredBase) where {F,S}
    n = length(d.pars)
    inbounds(x) = all(xj -> b.inbounds(xj), x)
    constℓ = n * b.constℓ
    varℓ() = n * b.varℓ()
    base = b.base ^ size(d.pars)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

# Same as PowerMeasure
@inline function logdensity(d::ProductMeasure{F,S,<:Fill}, x) where {F,S}
    d1 = d.f(first(d.pars))
    sum(xj -> logdensity(d1, xj), x)
end
