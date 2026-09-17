using Test

using MeasureBase
using MeasureBase: mspace_elsize, mspace_flatsize, NoMSpaceElementSize
using MeasureBase: preferred_stdmeasure, promote_stdmeasure, AnyStdMeasure, NoStdTransport
using MeasureBase: StdNormal, StdUniform, StdExponential, StdLogistic
using MeasureBase: Dirac, Lebesgue, Counting, LebesgueBase, CountingBase
using MeasureBase: mreshape, productmeasure, weightedmeasure, pushfwd, mbind, restrict
using IntervalSets: (..)
using StaticArrays: SVector, Size
using Static: static
using MeasureBase: size2length
using MeasureBase: setcartpower, ℝ, testvalue

_flat_iter(x::Number) = (x,)
_flat_iter(x::AbstractArray) = Iterators.flatten(map(_flat_iter, x))

struct _CustomStd <: MeasureBase.StdMeasure end

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
        @test @inferred(mspace_elsize(mreshape((StdNormal()^2)^6, (2, 3)))) == (2, 3)
        @test @inferred(mspace_flatsize(mreshape((StdNormal()^2)^6, (2, 3)))) isa NoMSpaceElementSize

        s = setcartpower(setcartpower(ℝ, 2), 3)
        @test @inferred(mspace_elsize(Lebesgue(s))) == (3,)
        @test @inferred(mspace_flatsize(Lebesgue(s))) == (2, 3)
        for μ in (StdNormal(), StdNormal()^3, (StdNormal()^2)^3, Dirac([1.0, 2.0]), Dirac(3.0))
            @test size2length(mspace_flatsize(μ)) == length(vec(collect(Iterators.flatten(_flat_iter(testvalue(μ))))))
        end

        @test @inferred(mspace_elsize(productmeasure((a = StdNormal(), b = StdUniform())))) isa NoMSpaceElementSize
        @test @inferred(mspace_flatsize(mbind(x -> StdNormal()^2, StdUniform()))) isa NoMSpaceElementSize
    end

    @testset "preferred_stdmeasure" begin
        for S in (StdNormal, StdUniform, StdExponential, StdLogistic)
            @test @inferred(preferred_stdmeasure(S())) === S
            @test @inferred(preferred_stdmeasure(S()^3)) === S
            @test @inferred(preferred_stdmeasure(weightedmeasure(0.1, S()))) === S
            @test @inferred(preferred_stdmeasure(pushfwd(exp, S()))) === S
            @test @inferred(preferred_stdmeasure(restrict(x -> x > 0, S()))) === S
        end

        @test @inferred(preferred_stdmeasure(Dirac(2.0))) === AnyStdMeasure
        @test @inferred(preferred_stdmeasure(Lebesgue())) <: NoStdTransport
        @test @inferred(preferred_stdmeasure(Counting())) <: NoStdTransport

        @test @inferred(preferred_stdmeasure(productmeasure((StdUniform(), StdNormal())))) === StdNormal
        @test @inferred(preferred_stdmeasure(productmeasure((a = StdUniform(), b = StdExponential())))) === StdExponential
        @test @inferred(preferred_stdmeasure(productmeasure((a = Dirac(1.0), b = StdUniform())))) === StdUniform
        @test @inferred(preferred_stdmeasure(productmeasure((a = Lebesgue(), b = StdUniform())))) <: NoStdTransport
        @test @inferred(preferred_stdmeasure(productmeasure(fill(StdLogistic(), 3)))) === StdLogistic
        @test @inferred(preferred_stdmeasure(productmeasure(()))) === AnyStdMeasure
        @test @inferred(preferred_stdmeasure(mbind(x -> StdNormal()^2, StdUniform()))) === StdUniform
    end

    @testset "promote_stdmeasure" begin
        @test @inferred(promote_stdmeasure(StdUniform, StdNormal)) === StdNormal
        @test @inferred(promote_stdmeasure(StdNormal, StdUniform)) === StdNormal
        @test @inferred(promote_stdmeasure(StdUniform, StdExponential)) === StdExponential
        @test @inferred(promote_stdmeasure(StdExponential, StdLogistic)) === StdLogistic
        @test @inferred(promote_stdmeasure(StdLogistic, StdNormal)) === StdNormal
        @test @inferred(promote_stdmeasure(StdLogistic, StdLogistic)) === StdLogistic
        @test @inferred(promote_stdmeasure(AnyStdMeasure, StdUniform)) === StdUniform
        @test @inferred(promote_stdmeasure(StdUniform, AnyStdMeasure)) === StdUniform
        @test @inferred(promote_stdmeasure(AnyStdMeasure, AnyStdMeasure)) === AnyStdMeasure
        @test @inferred(promote_stdmeasure(NoStdTransport{Int}, StdNormal)) === NoStdTransport{Int}
        @test @inferred(promote_stdmeasure(StdNormal, NoStdTransport{Int})) === NoStdTransport{Int}
        @test @inferred(promote_stdmeasure(NoStdTransport{Int}, AnyStdMeasure)) === NoStdTransport{Int}
        @test @inferred(promote_stdmeasure(AnyStdMeasure, NoStdTransport{Int})) === NoStdTransport{Int}
        @test @inferred(promote_stdmeasure(StdUniform, StdExponential, AnyStdMeasure, StdLogistic)) === StdLogistic
        @test @inferred(promote_stdmeasure(_CustomStd, StdUniform)) === StdUniform
        @test @inferred(promote_stdmeasure(_CustomStd, AnyStdMeasure)) === _CustomStd
        @test @inferred(preferred_stdmeasure(_CustomStd())) === _CustomStd
    end
end
