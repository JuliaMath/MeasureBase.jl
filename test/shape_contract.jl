using Test

using MeasureBase
using MeasureBase: mspace_elsize, mspace_flatsize, NoMSpaceElementSize
using MeasureBase: StdNormal, StdUniform, StdExponential, StdLogistic
using MeasureBase: Dirac, Lebesgue, Counting, LebesgueBase, CountingBase
using MeasureBase: mreshape, productmeasure, weightedmeasure, pushfwd, mbind, restrict
using IntervalSets: (..)
using StaticArrays: SVector, Size
using Static: static

@testset "shape contract" begin
    @testset "mspace_elsize and mspace_flatsize" begin
        for μ in (StdNormal(), StdUniform(), Lebesgue(), Lebesgue(0..1), Counting(), LebesgueBase(), CountingBase(), Dirac(1.5))
            @test @inferred(mspace_elsize(μ)) === ()
            @test @inferred(mspace_flatsize(μ)) === ()
        end

        @test @inferred(mspace_elsize(StdNormal()^3)) == (3,)
        @test @inferred(mspace_flatsize(StdNormal()^3)) == (3,)
        @test @inferred(mspace_elsize(StdNormal()^(2, 3))) == (2, 3)
        @test @inferred(mspace_flatsize(StdNormal()^(2, 3))) == (2, 3)

        @test @inferred(mspace_elsize((StdNormal()^3)^(2, 4))) == (2, 4)
        @test @inferred(mspace_flatsize((StdNormal()^3)^(2, 4))) == (3, 2, 4)
        @test @inferred(mspace_flatsize((StdNormal()^static(3))^static(2))) === Size(3, 2)
        @test @inferred(mspace_flatsize((StdNormal()^static(3))^2)) === (static(3), 2)

        @test @inferred(mspace_elsize(Dirac([1, 2]))) == (2,)
        @test @inferred(mspace_flatsize(Dirac([1, 2]))) == (2,)
        @test @inferred(mspace_flatsize(Dirac(SVector(1, 2)))) === Size(2)
        @test @inferred(mspace_elsize(Dirac([[1], [2]]))) == (2,)
        @test @inferred(mspace_flatsize(Dirac([[1], [2]])))  isa NoMSpaceElementSize
        @test @inferred(mspace_elsize(Dirac((a = 1, b = 2)))) isa NoMSpaceElementSize

        @test @inferred(mspace_elsize(weightedmeasure(0.3, StdNormal()^2))) == (2,)
        @test @inferred(mspace_flatsize(weightedmeasure(0.3, (StdNormal()^2)^3))) == (2, 3)
        @test @inferred(mspace_elsize(restrict(x -> x > 0, StdNormal()))) === ()
        @test @inferred(mspace_elsize(mreshape(StdNormal()^6, (2, 3)))) == (2, 3)
        @test @inferred(mspace_flatsize(mreshape(StdNormal()^6, (2, 3)))) == (2, 3)

        @test @inferred(mspace_elsize(productmeasure((a = StdNormal(), b = StdUniform())))) isa NoMSpaceElementSize
        @test @inferred(mspace_flatsize(mbind(x -> StdNormal()^2, StdUniform()))) isa NoMSpaceElementSize
    end
end
