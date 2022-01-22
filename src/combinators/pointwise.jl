export ⊙

struct PointwiseProductMeasure{M,L} <: AbstractMeasure
    measure::M
    likelihood::L

    function PointwiseProductMeasure(μ::M, ℓ::L) where {M,L}
        @assert static_hasmethod(logdensity_def, Tuple{L, gentype(μ)})
        return new{M,L}(μ, ℓ)
    end
end

function Base.show(io::IO, μ::PointwiseProductMeasure)
    io = IOContext(io, :compact => true)
    print(io, μ.measure, " ⊙ ", μ.likelihood)
end

function Base.show_unquoted(io::IO, μ::PointwiseProductMeasure, indent::Int, prec::Int)
    io = IOContext(io, :compact => true)
    if Base.operator_precedence(:*) ≤ prec
        print(io, "(")
        show(io, μ)
        print(io, ")")
    else
        show(io, μ)
    end
    return nothing
end

⊙(μ, ℓ) = pointwiseproduct(μ, ℓ)

@inline function logdensity_def(d::PointwiseProductMeasure, x)
    logdensity_def(d.measure, x) + logdensity_def(d.likelihood, x)
end

function gentype(d::PointwiseProductMeasure)
    @inbounds gentype(d.measure)
end

basemeasure(d::PointwiseProductMeasure) = @inbounds basemeasure(d.measure)
