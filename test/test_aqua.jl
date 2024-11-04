# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import MeasureBase

#Test.@testset "Package ambiguities" begin
#    Test.@test isempty(Test.detect_ambiguities(MeasureBase))
#end # testset

Test.@testset "Aqua tests" begin
    Aqua.test_all(MeasureBase, ambiguities = false)
end # testset
