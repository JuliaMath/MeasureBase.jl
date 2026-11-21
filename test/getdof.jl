using Test

using MeasureBase: getdof, check_dof, checked_arg
using MeasureBase: StdUniform, StdExponential, StdLogistic
using ChainRulesTestUtils: test_rrule
using Static: static

@testset "getdof" begin
    μ0 = StdExponential()
    x0 = rand(μ0)

    μ2 = StdExponential()^(2, 3)
    x2 = rand(μ2)

    @test @inferred(getdof(μ0)) === static(1)
    @test (check_dof(μ0, StdUniform()); true)
    @test_throws ArgumentError check_dof(μ2, μ0)
    test_rrule(check_dof, μ0, StdUniform())

    @test @inferred(checked_arg(μ0, x0)) === x0
    @test_throws ArgumentError checked_arg(μ0, x2)
    test_rrule(checked_arg, μ0, x0)

    @test @inferred(getdof(μ2)) == 6
    @test (check_dof(μ2, StdUniform()^(1, 6, 1)); true)
    @test_throws ArgumentError check_dof(μ2, μ0)
    test_rrule(check_dof, μ2, StdUniform()^(1, 6, 1))

    @test @inferred(checked_arg(μ2, x2)) === x2
    @test_throws ArgumentError checked_arg(μ2, x0)
    test_rrule(checked_arg, μ2, x2)

    @test @inferred(getdof((StdExponential()^3)^(static(0), static(0)))) === static(0)
end

@testset "fast_dof" begin
    using MeasureBase: fast_dof, NoFastDOF, NoDOF, StdNormal, Lebesgue, Dirac, Counting
    using MeasureBase: weightedmeasure, productmeasure, mbind, mcombine, transport_to, batched_rand_impl, GenContext
    using InverseFunctions: inverse

    # Static for primitive and standard measures and their wrappers:
    @test @inferred(fast_dof(StdNormal())) === static(1)
    @test @inferred(fast_dof(Lebesgue())) === static(1)
    @test @inferred(fast_dof(Dirac(0.3))) === static(0)
    @test @inferred(fast_dof(weightedmeasure(0.3, StdExponential()))) === static(1)
    @test fast_dof(Counting(Int)) isa NoDOF
    # Constant time for powers, products and combinations of such measures:
    @test @inferred(fast_dof(StdNormal()^static(3))) === static(3)
    @test fast_dof(StdNormal()^(2, 3)) == 6
    @test @inferred(fast_dof(productmeasure((a = StdNormal(), b = StdExponential()^static(2))))) === static(3)
    @test fast_dof(productmeasure([weightedmeasure(log(i), StdNormal()) for i in 1:4])) == 4
    @test fast_dof(mcombine(vcat, StdNormal()^2, StdUniform()^3)) == 5
    @test fast_dof(mbind(x -> StdNormal()^2, StdNormal(), vcat)) isa NoDOF

    # Products over arrays of marginals with several DOF each need a sum
    # over the marginals, transports and draws use the exact DOF:
    p = productmeasure([StdExponential()^2 for _ in 1:3])
    @test fast_dof(p) isa NoFastDOF
    @test getdof(p) == 6
    f = transport_to(StdNormal()^6, p)
    x = [[0.5, 1.0], [2.0, 0.1], [0.3, 0.7]]
    @test inverse(f)(f(x)) ≈ x
    X = batched_rand_impl(GenContext{Float64}(), p, (2,))
    @test size(X) == (2, 3, 2) && all(>(0), X)
    @test size(rand(p)) == (3,) && all(xi -> all(>(0), xi), rand(p))
end
