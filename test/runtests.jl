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
include("shape_contract.jl")
include("logdensities.jl")
include("structured_batches.jl")
include("batched_regressions.jl")
include("fixed_size_arrays.jl")
include("support_conventions.jl")
include("numtype.jl")
include("transport.jl")
include("transport_batched.jl")
include("smf.jl")
include("domains.jl")

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
include("combinators/product.jl")

include("rand.jl")
include("rand_batched.jl")

include("distributions/test_distributions.jl")

# Reactant only supports 64-bit Linux and macOS, and some of its
# dependencies break already during precompilation on other platforms,
# so it can't be a static test dependency:
if Sys.WORD_SIZE == 64 && (Sys.islinux() || Sys.isapple()) && isempty(VERSION.prerelease)
    import Pkg
    Base.identify_package("Reactant") === nothing && Pkg.add("Reactant")
    include("test_reactant.jl")
end

include("test_docs.jl")
