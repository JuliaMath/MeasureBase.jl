# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import MeasureBase

Test.@testset "Package ambiguities" begin
    Test.@test isempty(Test.detect_ambiguities(MeasureBase, recursive = true))
end # testset

Test.@testset "Aqua tests" begin
    Aqua.test_all(
        MeasureBase,
        # Only used by package extensions:
        stale_deps = (ignore = [:ArgCheck, :ArraysOfArrays],),
        # The check installs the package into a fresh environment, which
        # can't resolve while HeterogeneousComputing 0.3 is unregistered.
        # Enable it again once that is released:
        persistent_tasks = false,
    )
end # testset
