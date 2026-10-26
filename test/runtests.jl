using Test
using Base.Iterators: take
using Random
using LinearAlgebra
import LogarithmicNumbers

using MeasureBase
using MeasureBase.Interface: test_interface, test_smf

include("test_aqua.jl")

include("test_primitive.jl")
include("test_standard.jl")
include("test_basics.jl")

include("getdof.jl")
include("shape_contract.jl")
include("logdensities.jl")
include("structured_batches.jl")
include("parameter_batches.jl")
include("batched_regressions.jl")
include("static_variates.jl")
include("fixed_size_arrays.jl")
include("support_conventions.jl")
include("numtype.jl")
include("transport.jl")
include("transport_batched.jl")
include("smf.jl")
include("domains.jl")

include("measure_operators.jl")
include("public_api.jl")
include("contracts.jl")

include("rand.jl")
include("rand_batched.jl")
