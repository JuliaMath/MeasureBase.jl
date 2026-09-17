using Test

using MeasureBase
using MeasureBase: mspace_elsize, mspace_flatsize
using Distributions
using LinearAlgebra: I

@testset "shape contract for Distributions" begin
    @test @inferred(mspace_elsize(Normal(1, 2))) === ()
    @test @inferred(mspace_elsize(asmeasure(Normal(1, 2)))) === ()
    @test @inferred(mspace_flatsize(asmeasure(MvNormal(zeros(3), I(3))))) == (3,)
    @test @inferred(mspace_elsize(asmeasure(MvNormal(zeros(3), I(3)))^4)) == (4,)
    @test @inferred(mspace_flatsize(asmeasure(MvNormal(zeros(3), I(3)))^4)) == (3, 4)
end
