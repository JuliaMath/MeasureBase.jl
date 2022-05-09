struct ConditionalMeasure{M,C} <: AbstractMeasure
    parent::M
    constraint::C
end

"""
    (m::AbstractMeasure) | constraint

Return a new measure by constraining `m` to satisfy `constraint`.

Note that the form of `constraint` will vary depending on the structure of a
given measure. For example, a measure over `NamedTuple`s may allow `NamedTuple`
constraints, while another may require `constraint` to be a predicate or a
function returning a real number (in which case the constraint could be
considered as the zero-set of that function). 

At the time of this writing, invariants required of this function are not yet
settled. Specifically, there's the question of normalization. It's common for
conditional distributions to be normalized, but this can often not be expressed
in closed form, and can be very expensive to compute. For more general measures,
the notion of normalization may not even make sense.

Because of this, this interface is not yet stable, and users should expect
upcoming changes.
"""
Base.:|(μ::AbstractMeasure, constraint) = condition(μ, constraint)

condition(μ, constraint) = ConditionalMeasure(μ, constraint)

@inline basemeasure(cm::ConditionalMeasure) = basemeasure(cm.parent) | cm.constraint

# @generated function Base.:|(μ::ProductMeasure{NamedTuple{M,T}}, constraint::NamedTuple{N}) where {M,T,N}
#     newkeys = tuple(setdiff(M, N)...)
#     quote
#         mar = marginals(μ)
#         productmeasure(NamedTuple{$newkeys}(mar))
#     end
# end

function Base.:|(
    μ::ProductMeasure{NamedTuple{M,T}},
    constraint::NamedTuple{N},
) where {M,T,N}
    productmeasure(merge(marginals(μ), rmap(Dirac, constraint)))
end

function Pretty.tile(d::ConditionalMeasure)
    Pretty.pair_layout(Pretty.tile(d.parent), Pretty.tile(d.constraint), sep = " | ")
end
