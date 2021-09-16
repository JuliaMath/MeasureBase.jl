export ⊙

@concrete terse struct PointwiseProductMeasure{T} <: AbstractMeasure
    data :: T
end

Base.size(μ::PointwiseProductMeasure) = size(μ.data)

function Base.show(io::IO, μ::PointwiseProductMeasure)
    io = IOContext(io, :compact => true)
    print(io, join(string.(μ.data), " ⊙ "))
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

Base.length(m::PointwiseProductMeasure{T}) where {T} = length(m.data)

⊙(args...) = pointwiseproduct(args...)

@inline function logdensity(d::PointwiseProductMeasure, x)
    sum((logdensity(dⱼ, x) for dⱼ in d.data))
end

function sampletype(d::PointwiseProductMeasure) 
    @inbounds sampletype(first(d.data))
end

basemeasure(d::PointwiseProductMeasure) =  @inbounds basemeasure(first(d.data))
