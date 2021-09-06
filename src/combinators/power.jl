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

const PowerMeasure{F,T,N,A} = ProductMeasure{F,Fill{T,N,A}}

function Base.:^(μ::AbstractMeasure, dims::Integer...)
    return μ ^ dims
end

function Base.:^(μ::M, dims::NTuple{N,I}) where {M<:AbstractMeasure,N,I<:Integer}
    powermeasure(μ, dims)
end

# Same as PowerMeasure
function Base.show(io::IO, d::ProductMeasure{F,<:Fill}) where {F}
    io = IOContext(io, :compact => true)
    print(io, d.f(first(d.pars)), " ^ ", size(d.pars))
end

# Same as PowerMeasure{F,T,1} where {F,T}
function Base.show(io::IO, d::ProductMeasure{F,Fill{T,1,A}}) where {F,T,A}
    io = IOContext(io, :compact => true)
    print(io, d.f(first(d.pars)), " ^ ", size(d.pars)[1])
end

# sampletype(d::PowerMeasure{M,N}) where {M,N} = @inbounds Array{sampletype(first(marginals(d))), N}



params(d::ProductMeasure{F,<:Fill}) where {F} = params(first(marginals(d)))

params(::Type{P}) where {F,P<:ProductMeasure{F,<:Fill}} = params(D)

# basemeasure(μ::PowerMeasure) = @inbounds basemeasure(first(μ.data))^size(μ.data)

@inline basemeasure(d::PowerMeasure) = _basemeasure(d, (basemeasure(d.f(first(d.pars)))))

@inline _basemeasure(d::PowerMeasure, b) = b ^ size(d.pars)

@inline function _basemeasure(d::PowerMeasure, b::FactoredBase)
    n = length(d.pars)
    inbounds(x) = all(xj -> b.inbounds(xj), x)
    constℓ = n * b.constℓ
    varℓ() = n * b.varℓ()
    base = b.base ^ size(d.pars)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

function logdensity(d::PowerMeasure, x)
    d1 = d.f(first(d.pars))
    sum(xj -> logdensity(d1, xj), x)
end
