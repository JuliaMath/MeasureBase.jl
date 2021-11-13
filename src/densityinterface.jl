import DensityInterface

@inline DensityInterface.hasdensity(::AbstractMeasure) = true


function DensityInterface.logdensityof(m::AbstractMeasure, x)
    _logdensityof_impl(m, basemeasure(m), x, zero(Float64))
end

@inline function _logdensityof_impl(m::AbstractMeasure, β::AbstractMeasure, x, ℓ)
    Δℓ = logdensity(m, x)
    newℓ = ℓ + Δℓ
    m == β && return newℓ
    _logdensityof_impl(β, basemeasure(β), x, newℓ)
end


# MeasureBase.logpdf and Distributions.logpdf are currently implemented in MeasureTheory,
# so we'll leave this bit out for now:
#=
Base.@deprecate MeasureBase.logpdf(d::AbstractMeasure, x) DensityInterface.logdensityof(d, x)

function Distributions.logpdf(d::AbstractMeasure, x)
    mass(d) ≈ 1 || throw(ArgumentError("logpdf requires a normalized measure"))
    DensityInterface.logdensityof(d, x)
end

function Distributions.pdf(d::AbstractMeasure, x)
    mass(d) ≈ 1 || throw(ArgumentError("logpdf requires a normalized measure"))
    DensityInterface.logdensityof(d, x)
end
=#
