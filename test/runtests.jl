import Pkg
# Reactant, Enzyme and LLVM must come from the same resolved test environment.
const TEST_REACTANT = haskey(Pkg.project().dependencies, "Reactant") &&
    Sys.WORD_SIZE == 64 && (Sys.islinux() || Sys.isapple()) && isempty(VERSION.prerelease)

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
