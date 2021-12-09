struct ConditionalMeasure{M,C} <: AbstractMeasure
    parent::M 
    constraint::C
end

Base.:|(μ::AbstractMeasure, constraint) = ConditionalMeasure(μ, constraint)