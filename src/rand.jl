import Base

Base.rand(d::AbstractMeasure) = rand(Random.GLOBAL_RNG, Float64, d)

Base.rand(T::Type, μ::AbstractMeasure) = rand(Random.GLOBAL_RNG, T, μ)

Base.rand(rng::AbstractRNG, d::AbstractMeasure) = rand(rng, Float64, d)

@inline Random.rand!(d::AbstractMeasure, args...) = rand!(GLOBAL_RNG, d, args...)

# TODO: Make this work
# function Base.rand(rng::AbstractRNG, ::Type{T}, d::AbstractMeasure) where {T}
#     x = testvalue(d)
#     rand!(d, x)
# end

# struct ArraySlot{A,I}
#     arr::A
#     i::I
# end

# function rand!(rng::AbstractRNG, d::AbstractMeasure, x::ArraySlot)
#     x.arr[x.i...] = rand(rng, d)
# end
