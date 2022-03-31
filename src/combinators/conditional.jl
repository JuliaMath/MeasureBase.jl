struct ConditionalMeasure{M,C} <: AbstractMeasure
    parent::M 
    constraint::C
end

Base.:|(μ::AbstractMeasure, constraint) = ConditionalMeasure(μ, constraint)

@inline basemeasure(cm::ConditionalMeasure) = basemeasure(cm.parent) | cm.constraint