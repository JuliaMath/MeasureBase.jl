using Test

import MeasureBase

import Static
using Static: static
import FillArrays

@testset "static" begin
    @test 2 isa MeasureBase.IntegerLike
    @test static(2) isa MeasureBase.IntegerLike
    @test true isa MeasureBase.IntegerLike
    @test static(true) isa MeasureBase.IntegerLike

    @test @inferred(MeasureBase.one_to(7)) isa Base.OneTo
    @test @inferred(MeasureBase.one_to(7)) == 1:7
    @test @inferred(MeasureBase.one_to(static(7))) isa Static.SOneTo
    @test @inferred(MeasureBase.one_to(static(7))) == static(1):static(7)

    @test @inferred(MeasureBase.maybestatic_fill(4.2, (7,))) == FillArrays.Fill(4.2, 7)
    @test @inferred(MeasureBase.maybestatic_fill(4.2, (static(7),))) == FillArrays.Fill(4.2, 7)
    @test @inferred(MeasureBase.maybestatic_fill(4.2, (3, static(7)))) ==
          FillArrays.Fill(4.2, 3, 7)
    @test @inferred(MeasureBase.maybestatic_fill(4.2, (3:7,))) == FillArrays.Fill(4.2, (3:7,))
    @test @inferred(MeasureBase.maybestatic_fill(4.2, (static(3):static(7),))) ==
          FillArrays.Fill(4.2, (3:7,))
    @test @inferred(MeasureBase.maybestatic_fill(4.2, (3:7, static(2):static(5)))) ==
          FillArrays.Fill(4.2, (3:7, 2:5))

    @test MeasureBase.maybestatic_length(MeasureBase.one_to(7)) isa Int
    @test MeasureBase.maybestatic_length(MeasureBase.one_to(7)) == 7
    @test MeasureBase.maybestatic_length(MeasureBase.one_to(static(7))) isa Static.StaticInt
    @test MeasureBase.maybestatic_length(MeasureBase.one_to(static(7))) == static(7)
end
