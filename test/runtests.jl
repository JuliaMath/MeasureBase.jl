using Test
using Base.Iterators: take
using Random
using LinearAlgebra
import LogarithmicNumbers

using MeasureBase
using MeasureBase: test_interface, test_smf

include("test_aqua.jl")


include("test_primitive.jl")
include("test_standard.jl")
include("test_basics.jl")

include("getdof.jl")
include("transport.jl")
include("smf.jl")

include("test_mooncake.jl")

include("measure_operators.jl")

include("combinators/smart_constructors.jl")
include("combinators/weighted.jl")
include("combinators/superpose.jl")
include("combinators/transformedmeasure.jl")
include("combinators/reshape.jl")
include("combinators/implicitlymapped.jl")
include("combinators/combined.jl")
include("combinators/bind.jl")

include("distributions/test_distributions.jl")

include("test_docs.jl")
