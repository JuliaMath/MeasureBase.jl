# abstract type LogdensityComputation end

# struct LogpdfComputation end

# struct BasemeasureComputation{B} <: LogdensityComputation
#     base::B
# end

# struct LogdensityThunk{C,M,X,L}
#     computation::C
#     measure::M
#     x::X
#     ℓ::L
# end

# function reduce_step(thunk :: LogdensityThunk, callback=Returns(nothing))



# reduce_step(LogdensityThunk{M,B,X,L}, callback=Returns(nothing)) -> Union{LogdensityThunk, L}

# function logdensity_def(μ, b, x; cb=Returns(nothing))
#     f(thunk) = reduce_step(thunk, cb)
#     fix(f, LogdensityThunk(μ, b, x, zero(Float64)))
# end
