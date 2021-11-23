export CountingMeasure

struct CountingMeasure{X} <: PrimitiveMeasure end

function Base.show(io::IO, μ::CountingMeasure{X}) where {X}
    io = IOContext(io, :compact => true)
    print(io, "CountingMeasure(", X, ")")
end

CountingMeasure(X) = CountingMeasure{X}()

# gentype(::CountingMeasure{ℝ}) = Float64
# gentype(::CountingMeasure{ℝ₊}) = Float64
# gentype(::CountingMeasure{𝕀}) = Float64

gentype(::CountingMeasure) = Int

testvalue(μ::CountingMeasure{X}) where {X} = testvalue(X)

logdensity_def(::CountingMeasure, x) = zero(float(x))

# (::CountingMeaure)(s) = length(Set(s))
