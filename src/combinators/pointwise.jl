export ⊙

struct PointwiseProductMeasure{P,L} <: AbstractMeasure
    prior::P
    likelihood::L
end



iterate(p::PointwiseProductMeasure, i=1) = iterate((p.prior, p.likelihood), i)

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
    μ, ℓ = d
    logdensityof(ℓ.k(p), ℓ.x)
end

function gentype(d::PointwiseProductMeasure)
    gentype(d.prior)
end

@inbounds function insupport(d::PointwiseProductMeasure, p) 
    μ, ℓ = d
    insupport(μ, p) && insupport(ℓ.k(p), ℓ.x)
end

basemeasure(d::PointwiseProductMeasure, x) = d.prior

basemeasure(d::PointwiseProductMeasure) = basemeasure(d.prior)
