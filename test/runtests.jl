using Test
using Base.Iterators: take
using Random
using LinearAlgebra
import LogarithmicNumbers

using MeasureBase
using MeasureBase: test_interface, test_smf

include("test_aqua.jl")

include("static.jl")

include("test_basics.jl")

include("getdof.jl")
include("mspace.jl")
include("transport.jl")
include("smf.jl")

include("combinators/weighted.jl")
include("combinators/transformedmeasure.jl")
include("combinators/reshape.jl")
include("combinators/combined.jl")
include("combinators/bind.jl")

include("measure_operators.jl")

include("test_distributions.jl")

include("test_docs.jl")
