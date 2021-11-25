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

const PowerMeasure{M,N,R} = ProductMeasure{KernelReturns{M}, CartesianIndices{N, R}} 

Base.:^(μ::AbstractMeasure, ::Tuple{}) = μ

function Base.:^(μ::AbstractMeasure, dims::Integer...)
    return μ^dims
end

Base.show(io::IO, d::PowerMeasure) = print(io, d.f.value, " ^ ", size(d.pars))
Base.show(io::IO, d::PowerMeasure{M,1}) where {M} = print(io, d.f.value, " ^ ", length(d.pars))

# gentype(d::PowerMeasure{M,N}) where {M,N} = @inbounds Array{gentype(first(marginals(d))), N}

params(d::PowerMeasure) = params(first(marginals(d)))

params(::Type{P}) where {K,P<:PowerMeasure} = params(D)

# basemeasure(μ::PowerMeasure) = @inbounds basemeasure(first(μ.data))^size(μ.data)

# Same as PowerMeasure
@inline function basemeasure(d::PowerMeasure)
    _basemeasure(d, (basemeasure(d.f(first(d.pars)))))
end

# Same as PowerMeasure
@inline function _basemeasure(d::PowerMeasure, b)
    b ^ size(d.pars)
end

# Same as PowerMeasure
@inline function _basemeasure(d::PowerMeasure, b::FactoredBase)
    n = length(d.pars)
    inbounds(x) = all(b.inbounds, x)
    constℓ = n * b.constℓ
    varℓ() = n * b.varℓ()
    base = b.base^size(d.pars)
    FactoredBase(inbounds, constℓ, varℓ, base)
end

function basemeasure_depth(d::PowerMeasure) 
    return static(1) + basemeasure_depth(M)
end

function basemeasure_depth(::Type{P}) where {P<:PowerMeasure}
    return static(1) + basemeasure_depth(M)
end

# Same as PowerMeasure
@inline function logdensity_def(d::PowerMeasure, x)
    d1 = d.f(first(d.pars))
    sum(xj -> logdensity_def(d1, xj), x)
end
