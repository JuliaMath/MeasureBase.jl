using Test
using Base.Iterators: take
using Random
using LinearAlgebra
import LogarithmicNumbers

using MeasureBase
using MeasureBase: test_interface, test_smf

include("test_aqua.jl")

include("static.jl")

include("test_primitive.jl")
include("test_standard.jl")
include("test_basics.jl")

include("getdof.jl")
include("transport.jl")
include("smf.jl")

include("combinators/weighted.jl")
include("combinators/transformedmeasure.jl")
include("combinators/implicitlymapped.jl")

include("test_docs.jl")
