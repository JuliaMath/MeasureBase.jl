using Random
export ResettableRNG

export reset!
struct ResettableRNG{R,S} <: Random.AbstractRNG
    rng::R
    seed::S
end

function Base.show(io::IO, r::ResettableRNG)
    io = IOContext(io, :compact => true)
    print(io, "ResettableRNG(::", constructor(r.rng), ", ", r.seed, ")")
end

function reset!(r::ResettableRNG)
    @info "Calling reset!"
    Random.seed!(r.rng, r.seed)
end


# for f in [
#     :(Base.rand)
#     :(Base.randn)
#     :(Random.rand!)
#     :(Random.randcycle)
#     :(Random.randcycle!)
#     :(Random.randexp)
#     :(Random.randexp!)
#     :(Random.randn!)
#     :(Random.randperm)
#     :(Random.randperm!)
#     :(Random.randstring)
#     :(Random.randsubseq)
#     :(Random.randsubseq!)
#     :(Random.shuffle)
#     :(Random.shuffle!)
#     :(Random.seed!)
# ]    
#     @eval $f(r::ResettableRNG, args...) = $f(r.rng, args...)
# end


# Base.rand(r::ResettableRNG, d::AbstractMeasure) = rand(r.rng, d)
# Base.rand(r::ResettableRNG, ::Type{T}, d::AbstractMeasure) where {T} = rand(r.rng, T, d)
# Base.rand(r::ResettableRNG) = rand(r.rng, Float64)

import Random
using Random: randexp

for T in isbits_subtypes(Real)
    @eval begin
        
        function Base.rand(r::ResettableRNG, ::Type{$T}) 
            rand(r.rng, $T)
        end

        function Base.randn(r::ResettableRNG, ::Type{$T}) 
            randn(r.rng, $T)
        end

        function Random.randexp(r::ResettableRNG, ::Type{$T})
            randexp(r.rng, $T)
        end
    end
end

function Base.iterate(r::ResettableRNG)
    r = deepcopy(r)
    reset!(r)
    return (rand(r), nothing)
end

Base.iterate(r::ResettableRNG, _) = (rand(r), nothing)
Base.IteratorSize(r::ResettableRNG) = Base.IsInfinite()

function Random.Sampler(r::Type{R}, s::Random.Sampler, rep::Random.Repetition) where {R<:ResettableRNG}
    return Random.Sampler(r.rng, s, r)
end

function Base.rand(r::ResettableRNG, sp::Random.Sampler)
    rand(r.rng, sp)
end
