using Test

using MeasureBase: getdof, check_dof, checked_arg
using MeasureBase: StdUniform, StdExponential, StdLogistic
using ChainRulesTestUtils: test_rrule
using Static: static


@testset "getdof" begin
    μ0 = StdExponential()
    x0 = rand(μ0)

    μ2 = StdExponential()^(2,3)
    x2 = rand(μ2)

    @test @inferred(getdof(μ0)) === static(1)
    @test (check_dof(μ0, StdUniform()); true)
    @test_throws ArgumentError check_dof(μ2, μ0)
    test_rrule(check_dof, μ0, StdUniform())

    @test @inferred(checked_arg(μ0, x0)) === x0
    @test_throws ArgumentError checked_arg(μ0, x2)
    test_rrule(checked_arg, μ0, x0)

    @test @inferred(getdof(μ2)) == 6
    @test (check_dof(μ2, StdUniform()^(1,6,1)); true)
    @test_throws ArgumentError check_dof(μ2, μ0)
    test_rrule(check_dof, μ2, StdUniform()^(1,6,1))

    @test @inferred(checked_arg(μ2, x2)) === x2
    @test_throws ArgumentError checked_arg(μ2, x0)
    test_rrule(checked_arg, μ2, x2)

    @test @inferred(getdof((StdExponential()^3)^(static(0),static(0)))) === static(0)
end
