using Test

using MeasureBase
using MeasureBase: mreshape, Reshape
using InverseFunctions: inverse
using ChangesOfVariables: with_logabsdet_jacobian

@testset "reshape" begin
    μ = StdNormal()^6
    ν = mreshape(μ, 2, 3)
    X = reshape(collect(0.1:0.1:0.6), 2, 3)

    @test size(rand(ν)) == (2, 3)
    @test logdensityof(ν, X) ≈ logdensityof(μ, vec(X))
    Xs = rand(2, 3, 4)
    @test logdensities(ν, Xs) ≈ [logdensityof(μ, vec(Xs[:, :, i])) for i in 1:4]
    @test logdensityof(ν, fill(Inf, 2, 3)) == -Inf
    @test logdensityof(mreshape(StdUniform()^4, (2, 2)), [0.5 1.5; 0.2 0.3]) == -Inf
    @test_throws ArgumentError logdensityof(ν, zeros(3, 2))

    f = transport_to(StdUniform()^6, ν)
    @test f(X) ≈ transport_to(StdUniform()^6, μ)(vec(X))
    @test inverse(f)(f(X)) ≈ X

    g = Reshape((2, 3), (6,))
    @test g(vec(X)) == X
    @test inverse(g)(X) == vec(X)
    @test with_logabsdet_jacobian(g, vec(X)) == (X, 0.0)
end
