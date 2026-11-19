# Resolve Reactant before loading Enzyme or LLVM. Adding it can change
# their versions. Reactant is only tested on supported stable platforms.
const TEST_REACTANT = Sys.WORD_SIZE == 64 && (Sys.islinux() || Sys.isapple()) && isempty(VERSION.prerelease)
if TEST_REACTANT && Base.identify_package("Reactant") === nothing
    import Pkg
    Pkg.add("Reactant")
end

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
include("root_relations.jl")
include("numtype.jl")
include("transport.jl")
include("transport_batched.jl")
include("smf.jl")
include("moments.jl")
include("domains.jl")

include("test_mooncake.jl")
include("test_zygote.jl")

include("measure_operators.jl")
include("public_api.jl")
include("contracts.jl")
include("extension_recipes.jl")

include("combinators/smart_constructors.jl")
include("combinators/weighted.jl")
include("combinators/normalized.jl")
include("combinators/truncated.jl")
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

TEST_REACTANT && include("test_reactant.jl")

@testset "Loaded extension ambiguities" begin
    import Pkg
    names = keys(Pkg.TOML.parsefile(pkgdir(MeasureBase, "Project.toml"))["extensions"])
    extensions = Module[]
    for name in names
        name == "MeasureBaseReactantExt" && !TEST_REACTANT && continue
        ext = Base.get_extension(MeasureBase, Symbol(name))
        @test ext isa Module
        isnothing(ext) || push!(extensions, ext)
    end
    @test isempty(Test.detect_ambiguities(MeasureBase, extensions...; recursive = true))
end

include("test_docs.jl")
