using Test

using MeasureBase
using MeasureBase: mspace_elsize, mspace_flatsize, preferred_stdmeasure, NoStdTransport
using MeasureBase: StdNormal, StdUniform, StdExponential, StdLogistic, Dirac, productmeasure
using Distributions
using LinearAlgebra: I

@testset "shape contract for Distributions" begin
    @test @inferred(mspace_elsize(Normal(1, 2))) === ()
    @test @inferred(mspace_elsize(asmeasure(Normal(1, 2)))) === ()
    @test @inferred(mspace_flatsize(asmeasure(MvNormal(zeros(3), I(3))))) == (3,)
    @test @inferred(mspace_elsize(asmeasure(MvNormal(zeros(3), I(3)))^4)) == (4,)
    @test @inferred(mspace_flatsize(asmeasure(MvNormal(zeros(3), I(3)))^4)) == (3, 4)

    @test @inferred(preferred_stdmeasure(Normal(1, 2))) === StdNormal
    @test @inferred(preferred_stdmeasure(asmeasure(Normal(1, 2)))) === StdNormal
    @test @inferred(preferred_stdmeasure(3 + 2 * Normal())) === StdNormal
    @test @inferred(preferred_stdmeasure(Uniform(1, 2))) === StdUniform
    @test @inferred(preferred_stdmeasure(Exponential(2.0))) === StdExponential
    @test @inferred(preferred_stdmeasure(Logistic(1, 2))) === StdLogistic
    @test @inferred(preferred_stdmeasure(Beta(2, 3))) === StdUniform
    @test @inferred(preferred_stdmeasure(truncated(Normal(), 0, 1))) === StdUniform
    @test @inferred(preferred_stdmeasure(MvNormal(zeros(2), I(2)))) === StdNormal
    @test @inferred(preferred_stdmeasure(Dirichlet([1.0, 2.0]))) === StdUniform
    @test @inferred(preferred_stdmeasure(Poisson(3))) <: NoStdTransport
    @test @inferred(preferred_stdmeasure(StandardDist{Normal}(3))) === StdNormal
    @test @inferred(preferred_stdmeasure(StandardDist{Uniform}())) === StdUniform

    @test @inferred(preferred_stdmeasure(productmeasure((asmeasure(Beta(2, 3)), asmeasure(Normal()))))) === StdNormal
    @test @inferred(preferred_stdmeasure(productmeasure((a = Dirac(1.0), b = asmeasure(Beta(2, 3)))))) === StdUniform
    @test @inferred(preferred_stdmeasure(productmeasure((a = asmeasure(Poisson(2)), b = asmeasure(Beta(2, 3)))))) <: NoStdTransport
    @test @inferred(preferred_stdmeasure(productmeasure([asmeasure(Normal(i, 1)) for i in 1:3]))) === StdNormal
end
