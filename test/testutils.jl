# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Allocations of `f(args...)` after a warm-up call, measured inside a
# function since `@allocated` at top level reports a boxed result on
# Julia 1.10:
function allocations_of(f::F, args::Vararg{Any,N}) where {F,N}
    f(args...)
    @allocated f(args...)
end
