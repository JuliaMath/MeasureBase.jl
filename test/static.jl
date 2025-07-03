using Test

import MeasureBase

import Static
using StaticArrays
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

    @test @inferred(MeasureBase.fill_with(4.2, (7,))) == FillArrays.Fill(4.2, 7)
    @test @inferred(MeasureBase.fill_with(4.2, (static(7),))) == FillArrays.Fill(4.2, 7)
    @test @inferred(MeasureBase.fill_with(4.2, (3, static(7)))) ==
          FillArrays.Fill(4.2, 3, 7)
    @test @inferred(MeasureBase.fill_with(4.2, (3:7,))) == FillArrays.Fill(4.2, (3:7,))
    @test @inferred(MeasureBase.fill_with(4.2, (static(3):static(7),))) ==
          FillArrays.Fill(4.2, (3:7,))
    @test @inferred(MeasureBase.fill_with(4.2, (3:7, static(2):static(5)))) ==
          FillArrays.Fill(4.2, (3:7, 2:5))

    @test MeasureBase.maybestatic_length(MeasureBase.one_to(7)) isa Int
    @test MeasureBase.maybestatic_length(MeasureBase.one_to(7)) == 7
    @test MeasureBase.maybestatic_length(MeasureBase.one_to(static(7))) isa Static.StaticInt
    @test MeasureBase.maybestatic_length(MeasureBase.one_to(static(7))) == static(7)
end

@testset "maybestatic_size" begin
    # Test regular array
    arr = rand(MersenneTwister(123), 3, 4)
    @test MeasureBase.maybestatic_size(arr) == (3, 4)
    
    # Test static array
    static_arr = SMatrix{2,2}([1 2; 3 4])
    @test MeasureBase.maybestatic_size(static_arr) == (static(2), static(2))
    
    # Test mixed static/dynamic array
    mixed = zeros(static(2), 3)  # Create a matrix with static first dimension
    size_result = MeasureBase.maybestatic_size(mixed)
    @test size_result[1] isa Static.StaticInt
    @test size_result[2] isa Int
    @test size_result == (static(2), 3)
end