using Test
using MeasureBase
using InverseFunctions: inverse

# Checks the batching, support and precision contracts of measure `μ` on
# the variates `xs`, `x_out` (if given) must lie outside of the support:
function test_measure_contracts(μ, xs::AbstractVector, x_out; transport::Bool = true)
    @testset "contracts of $μ" begin
        _f32(x::Real) = Float32(x)
        _f32(x::AbstractArray) = Float32.(x)
        _f32(x::Union{Tuple,NamedTuple}) = map(_f32, x)
        xs32 = map(_f32, xs)

        ld = map(x -> logdensityof(μ, x), xs)
        @test all(isfinite, ld)
        @test logdensities(μ, xs) ≈ ld
        @test isempty(logdensities(μ, xs[1:0]))
        @test logdensityof(μ, first(xs32)) isa Float32
        @test eltype(logdensities(μ, xs32)) == Float32
        if !isnothing(x_out)
            @test logdensityof(μ, x_out) == -Inf
            @test logdensities(μ, [x_out, first(xs)]) ≈ [-Inf, first(ld)]
        end

        if transport
            f = transport_to(StdUniform, μ)
            ys = map(f, xs)
            @test all(map((a, b) -> all(a .≈ b), f.(xs), ys))
            @test all(map((a, b) -> all(map(≈, a, b)), inverse(f).(ys), xs))
            @test isempty(f.(xs[1:0]))
            @test eltype(f(first(xs32))) == Float32
            isnothing(x_out) || @test any(isnan, f(x_out))
        end
    end
end

test_measure_contracts(μ, x_out; kwargs...) =
    test_measure_contracts(μ, [rand(μ) for _ in 1:5], x_out; kwargs...)

@testset "measure contracts" begin
    test_measure_contracts(StdNormal(), nothing)
    test_measure_contracts(StdLogistic(), nothing)
    test_measure_contracts(StdExponential(), -1.0)
    test_measure_contracts(StdUniform(), 1.5)
    test_measure_contracts(StdExponential()^3, [1.0, -1.0, 1.0])
    test_measure_contracts(StdExponential()^(2, 2), [1.0 -1.0; 1.0 1.0])
    test_measure_contracts(weightedmeasure(0.5, StdExponential()), -1.0)
    test_measure_contracts(half(StdNormal()), -1.0)
    test_measure_contracts(pushfwd(exp, StdNormal()), nothing)
    test_measure_contracts(productmeasure((StdNormal(), StdExponential())), (0.5, -1.0))
    test_measure_contracts(productmeasure((a = StdNormal(), b = StdUniform())), (a = 0.5, b = 2.0))
    test_measure_contracts(mcombine(vcat, StdNormal()^2, StdExponential()^2), [0.1, 0.2, -1.0, 1.0])
    test_measure_contracts(Dirac(2.0), 3.0; transport = false)
    test_measure_contracts(superpose(StdUniform(), StdExponential()), -1.0; transport = false)
    test_measure_contracts(restrict(>(0), StdNormal()), [0.5, 1.0, 2.0], -1.0; transport = false)
    test_measure_contracts(mintegrate_exp(x -> -x, StdExponential()), [0.5, 1.0, 2.0], -1.0; transport = false)

    # Known gaps:
    @test_broken logdensityof(pushfwd(exp, StdNormal()), -1.0) == -Inf
    let P = productmeasure([StdNormal(), StdExponential()])
        @test_broken isempty(logdensities(P, rand(P^0)))
    end
end
