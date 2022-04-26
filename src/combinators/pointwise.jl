export ⊙

struct PointwiseProductMeasure{P,L} <: AbstractMeasure
    prior::P
    likelihood::L
end

function Base.show(io::IO, μ::PointwiseProductMeasure)
    io = IOContext(io, :compact => true)
    print(io, μ.prior, " ⊙ ", μ.likelihood)
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

@inline function logdensity_def(d::PointwiseProductMeasure, p)
    (μ, ℓ) = (d.prior, d.likelihood)
    logdensity_def(μ, p) + logdensityof(ℓ.k(p), ℓ.x)
end

function gentype(d::PointwiseProductMeasure)
    gentype(d.prior)
end

basemeasure(d::PointwiseProductMeasure, x) = basemeasure(d.prior, x)

basemeasure(d::PointwiseProductMeasure) = basemeasure(d.prior)

insupport(d::PointwiseProductMeasure, x) = insupport(d.prior, x)