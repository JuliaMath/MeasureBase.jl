using Test
using Base.Iterators: take
using Random
using LinearAlgebra
import LogarithmicNumbers

using MeasureBase
using MeasureBase: test_interface, test_smf

using Aqua
@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(MeasureBase, ambiguities = false)
    # Aqua.test_ambiguities(MeasureBase)
end

using JET
# @testset "Code linting (JET.jl)" begin
#     JET.test_package(MeasureBase; target_defined_modules = true)
# end

# include("test_aqua.jl")

include("static.jl")
include("domains.jl")

include("test_primitive.jl")
include("test_standard.jl")
include("test_basics.jl")

include("getdof.jl")
include("transport.jl")
include("smf.jl")

include("combinators/weighted.jl")
include("combinators/transformedmeasure.jl")
include("combinators/implicitlymapped.jl")
include("combinators/conditional.jl")
include("combinators/half.jl")
include("test_docs.jl")
