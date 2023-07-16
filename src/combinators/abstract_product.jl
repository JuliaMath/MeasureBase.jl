"""
    marginals(μ::AbstractMeasure)

Returns the marginals measures of `μ` as a collection of measures.

The kind of marginalization implied by `marginals` depends on the
type of `μ`.

`μ` may be a power of a measure or a product of measures, but other
types of measures may support `marginals` as well.
"""
function marginals end
export marginals


"""
    abstract type AbstractProductMeasure

Abstact type for products of measures.

[`marginals(μ::AbstractProductMeasure)`](@ref) returns the collection of
measures that `μ` is the product of.
"""
abstract type AbstractProductMeasure <: AbstractMeasure end
export AbstractProductMeasure

function Pretty.tile(μ::AbstractProductMeasure)
    result = Pretty.literal("ProductMeasure(")
    result *= Pretty.tile(marginals(μ))
    result *= Pretty.literal(")")
end

massof(m::AbstractProductMeasure) = prod(massof, marginals(m))

Base.:(==)(a::AbstractProductMeasure, b::AbstractProductMeasure) = marginals(a) == marginals(b)
Base.isapprox(a::AbstractProductMeasure, b::AbstractProductMeasure; kwargs...) = isapprox(marginals(a), marginals(b); kwargs...)


# # ToDo: Do we want this? It's not so clear what the semantics of `length` and `size`
# # for measures should be, in general:
# Base.length(μ::AbstractProductMeasure) = length(marginals(μ))
# Base.size(μ::AbstractProductMeasure) = size(marginals(μ))
