using Test

using MeasureBase
using MeasureBase: pushfwd, StdUniform, StdExponential, StdLogistic
using MeasureBase: pushfwd, PushforwardMeasure
using MeasureBase: transport_to, unsafe_logdensityof
using MeasureBase: productmeasure, mbind, superpose, SpikeMixture
import AffineMaps
import Zygote
using InverseFunctions: inverse
using ChangesOfVariables: with_logabsdet_jacobian
import Statistics: var
using DensityInterface: logdensityof
using LogExpFunctions
using SpecialFunctions: erfc, erfcinv
import InverseFunctions: inverse, FunctionWithInverse, setinverse
using IrrationalConstants: invsqrt2, sqrt2
using IntervalSets
import ChangesOfVariables: with_logabsdet_jacobian
using MeasureBase.Interface: transport_to, test_transport

Φ(z) = erfc(-z * invsqrt2) / 2
Φinv(p) = -erfcinv(2 * p) * sqrt2

with_logabsdet_jacobian(f::FunctionWithInverse, x) = with_logabsdet_jacobian(f.f, x)

with_logabsdet_jacobian(::typeof(Φ), z) = (Φ(z), logdensityof(StdNormal(), z))

function with_logabsdet_jacobian(::typeof(Φinv), p)
    z = Φinv(p)
    (z, -logdensityof(StdNormal(), z))
end

inverse(::typeof(Φ)) = Φinv
inverse(::typeof(Φinv)) = Φ

var(::StdNormal) = 1.0
var(::StdExponential) = 1.0
var(::StdUniform) = 1 / 12
var(::StdLogistic) = π^2 / 3

function test_pushfwd(f, μ, ν_ref)
    @testset "pushfwd($f, $μ)" begin
        @inferred(pushfwd(f, μ))
        ν = pushfwd(f, μ)

        test_transport(ν, ν_ref)
        test_transport(ν_ref, ν)

        y = rand(ν_ref)
        @test isapprox(@inferred(logdensityof(ν, y)), logdensityof(ν_ref, y), atol = 1e-10)
        @test isapprox(
            @inferred(unsafe_logdensityof(ν, y)),
            unsafe_logdensityof(ν_ref, y),
            atol = 1e-10,
        )

        @test isapprox(var(rand(ν^(10^5))), var(ν_ref), rtol = 0.05)

        @test abs(transport_to(StdLogistic(), ν)(y)) ≈
              abs(transport_to(StdLogistic(), ν_ref)(y))
    end
end

@testset "transformedmeasure.jl" begin
    # (f, μ, ν_ref), so that pushfwd(f, μ) ≅ ν_ref
    triples = [
        ((-) ∘ log, StdUniform(), StdExponential())
        (exp ∘ (-), StdExponential(), StdUniform())
        (logit, StdUniform(), StdLogistic())
        (logistic, StdLogistic(), StdUniform())
        (Φ, StdNormal(), StdUniform())
        (Φinv, StdUniform(), StdNormal())
    ]

    for (f, μ, ν_ref) in triples
        test_pushfwd(f, μ, ν_ref)
        test_pushfwd(setinverse(f, inverse(f)), μ, ν_ref)
    end

    @testset "Pushforward-of-pushforward" begin
        for (f, μ, ν_ref) in triples
            finv = inverse(f)
            test_pushfwd(finv, pushfwd(f, μ), μ)
        end
    end
end

using IrrationalConstants: loghalf

@testset "Half" begin
    μ = Half(StdNormal())

    # Test logdensityof for positive values
    @test logdensityof(μ, 1.0) ≈ logdensityof(StdNormal(), 1.0) - loghalf

    # Test logdensityof for negative values (should return -Inf)
    @test logdensityof(μ, -1.0) == -Inf

    # Test logdensityof for zero
    @test logdensityof(μ, 0.0) ≈ logdensityof(StdNormal(), 0.0) - loghalf

    @test logdensities(μ, [-1.0, 1.0]) ≈ [-Inf, logdensityof(μ, 1.0)]
    @test MeasureBase.insupportof(μ, -1.0) == false && MeasureBase.insupportof(μ, 1.0) == true
    @test logdensityof(Half(2 * StdNormal()), 1.0) ≈ log(2) + logdensityof(μ, 1.0)

    # Half requires scalar variates and parents symmetric about zero:
    @test_throws ArgumentError Half(StdExponential())
    @test_throws ArgumentError Half(StdUniform())
    @test_throws ArgumentError Half(StdNormal()^3)

    # The fold would double the mass of atoms at zero, parents with atoms
    # are refused, continuous mixtures fold:
    @test_throws ArgumentError Half(Dirac(0.0))
    @test_throws ArgumentError Half(SpikeMixture(StdNormal(), 0.5))
    @test_throws ArgumentError Half(superpose(StdNormal(), Dirac(0.0)))
    φ(x) = exp(logdensityof(StdNormal(), x))
    mix = superpose(0.5 * pushfwd(AffineMaps.Add(-1.0), StdNormal()), 0.5 * pushfwd(AffineMaps.Add(1.0), StdNormal()))
    @test logdensityof(Half(mix), 0.7) ≈ log(φ(1.7) + φ(-0.3))

    # Half of a weighted parent has its mass, its transports map the
    # normalized law (Φ⁻¹(3/4) and 2Φ(1/2) - 1), also within truncations:
    h = Half(2 * StdNormal())
    @test massof(h) ≈ 2 && smf(h, 0.5) ≈ 2 * (2 * Φ(0.5) - 1) && smf(h, -1.0) == 0
    @test transport_to(h, StdUniform())(0.5) ≈ Φinv(0.75)
    @test transport_to(StdUniform(), h)(0.5) ≈ 2 * Φ(0.5) - 1
    @test transport_to(h, StdUniform())(0.5f0) isa Float32
    test_transport(h, StdUniform())
    for (a, b) in ((0.2, 1.5), (3.0, Inf))
        t, t1 = truncate(h, a .. b), truncate(Half(StdNormal()), a .. b)
        @test massof(t) ≈ 2 * massof(t1)
        @test transport_to(t, StdUniform())(0.3) ≈ transport_to(t1, StdUniform())(0.3)
        @test transport_to(StdUniform(), t)(a + 0.1) ≈ transport_to(StdUniform(), t1)(a + 0.1)
    end
end

using Test
using MeasureBase
using MeasureBase: gettransform
using InverseFunctions
using StaticArrays
using ChangesOfVariables
using ChangesOfVariables: setladj

@testset "PushforwardMeasure" begin
    # Test basic pushforward construction
    μ = StdNormal()
    f = exp
    ν = pushfwd(f, μ)

    @test ν isa PushforwardMeasure
    @test MeasureBase.paramsof(ν) == NamedTuple() && MeasureBase.paramnames(ν) == ()

    # TODO: How do we access the parent?
    # @test parent(ν) === μ
    @test gettransform(ν) === f

    # Test logdensity with AdaptRootMeasure
    x = 0.5
    y = f(x)
    ld = logdensityof(ν, y)
    @test !isnan(ld)
    @test isfinite(ld)

    # Test logdensity with PushfwdRootMeasure
    ν_no_corr = pushfwd(f, μ, PushfwdRootMeasure())
    ld_no_corr = logdensityof(ν_no_corr, y)
    @test !isnan(ld_no_corr)
    @test isfinite(ld_no_corr)

    # Test non-bijective pushforward
    g = x -> x^2
    @test MeasureBase.getdof(pushfwd(g, μ)) isa MeasureBase.NoDOF
    @test_throws ArgumentError logdensityof(pushfwd(g, μ), 1.0)
    @test_throws ArgumentError logdensityof(pushfwd(g, μ, PushfwdRootMeasure()), 1.0)

    # Test edge cases for _combine_logd_with_ladj
    @test MeasureBase._combine_logd_with_ladj(-Inf, Inf) == -Inf  # Zero density wins
    @test MeasureBase._combine_logd_with_ladj(1.0, -Inf) == -Inf
    @test logdensityof(pushfwd(log, StdExponential()), -Inf) == -Inf

    # Test composition of pushforwards
    h = x -> exp(x)
    ν_comp = pushfwd(h, ν)
    @test ν_comp isa PushforwardMeasure

    # TODO: How do we access the parent?
    # @test parent(ν_comp) === μ

    # Test identity pushforward
    ν_id = pushfwd(identity, μ)
    @test ν_id === μ

    # Test rootmeasure behavior
    @test rootmeasure(ν) === rootmeasure(μ)  # AdaptRootMeasure
    @test rootmeasure(ν_no_corr) isa PushforwardMeasure  # PushfwdRootMeasure

    # Pushforwards adapting the root have the origin's root as their base,
    # pushforwards of the root pull weights out of the pushed-forward base:
    @test basemeasure(ν) === rootmeasure(μ)
    @test basemeasure(ν_no_corr) isa WeightedMeasure{<:Any,<:PushforwardMeasure}
    @test basemeasure(ν_no_corr).base.style isa PushfwdRootMeasure

    # Both styles describe the same measure:
    @test logdensity_rel(ν, ν_no_corr, 2.0) ≈ 0 atol = 1e-12
    @test logdensity_rel(ν_no_corr, ν, 2.0) ≈ 0 atol = 1e-12
    @test logdensity_rel(ν, rootmeasure(μ), 2.0) ≈ logdensityof(μ, log(2.0)) - log(2.0)
    E = Base.BroadcastFunction(exp)
    νv, νv_no_corr = pushfwd(E, μ^2), pushfwd(E, μ^2, PushfwdRootMeasure())
    @test logdensity_rel(νv, νv_no_corr, [2.0, 3.0]) ≈ 0 atol = 1e-12
    @test logdensities(νv, [2.0 1.0; 3.0 0.5]) ≈ [logdensityof(νv, [2.0, 3.0]), logdensityof(νv, [1.0, 0.5])]

    # Rank-changing pushforwards reshape the root:
    mr = MeasureBase.mreshape(μ^6, (2, 3))
    @test basemeasure(mr) == MeasureBase.LebesgueBase()^(2, 3)
    xr = reshape(collect(0.1:0.1:0.6), 2, 3)
    @test logdensityof(mr, xr) ≈ logdensityof(μ^6, vec(xr))

    # Mixtures of pushforward location families:
    shifted = pushfwd(setinverse(setladj(x -> x + 1, zero), setladj(y -> y - 1, zero)), μ)
    s = superpose([μ, shifted])
    @test logdensityof(s, 0.3) ≈ log(exp(logdensityof(μ, 0.3)) + exp(logdensityof(μ, -0.7)))
    @test sum(x -> exp(logdensityof(s, x)), -12:0.01:13) * 0.01 ≈ 2

    # Test massof
    # TODO: mass interface is very incomplete
    # @test massof(ν) == massof(μ)

    # Test rand
    @test rand(ν) isa Real
    @test insupportof(ν, rand(ν)) != false

    # Test pullback
    pb = pullbck(f, ν)
    @test pb isa PushforwardMeasure
    @test logdensityof(pb, y) ≈ logdensityof(μ, y)

    # Test deprecated pullback
    @test_deprecated pullback(f, μ)
end

@testset "pushforwards of counting measures" begin
    # Injective maps keep the masses of atoms, without a volume element:
    D = MeasureBase.Dirac
    a = superpose(2D(1.0), 3D(2.0))
    p = pushfwd(AffineMaps.Mul(2.0), a)
    @test p.style isa PushfwdRootMeasure
    @test logdensityof(p, 2.0) ≈ log(2) && logdensityof(p, 4.0) ≈ log(3)
    @test logdensityof(p, 3.0) == -Inf
    @test sum(y -> exp(logdensityof(p, y)), 0.0:0.5:6.0) ≈ massof(a)
    @test pushfwd(AffineMaps.Mul(2.0), Counting(ℤ)).style isa PushfwdRootMeasure
    @test pushfwd(AffineMaps.Mul(2.0), StdNormal()).style isa AdaptRootMeasure
end

@testset "PushFwdStyle types" begin
    @test AdaptRootMeasure() isa PushFwdStyle
    @test PushfwdRootMeasure() isa PushFwdStyle
    @test MeasureBase.WithVolCorr === AdaptRootMeasure
    @test MeasureBase.NoVolCorr === PushfwdRootMeasure

    @testset "output size of pushforwards of tuple products" begin
        Pt = productmeasure((StdNormal(), StdUniform()^2))
        ν = pushfwd(x -> vcat(x[1], x[2]), Pt)
        @test MeasureBase.mspace_flatsize(ν) == (3,)
        @test MeasureBase.mspace_ndims(ν) == 1
        @test MeasureBase.mspace_flatsize(pushfwd(x -> x, mbind(x -> StdNormal()^(x > 0 ? 1 : 2), StdNormal()))) isa MeasureBase.NoMSpaceElementSize
    end

    @testset "construction inside differentiated functions" begin
        flat(x::Tuple) = vcat(x[1], x[2])
        unflat(v::AbstractVector) = (v[1], v[2:end])
        InverseFunctions.inverse(::typeof(flat)) = unflat
        InverseFunctions.inverse(::typeof(unflat)) = flat
        ChangesOfVariables.with_logabsdet_jacobian(::typeof(flat), x) = (flat(x), zero(eltype(x[2])))
        ChangesOfVariables.with_logabsdet_jacobian(::typeof(unflat), v) = (unflat(v), zero(eltype(v)))
        Pt = productmeasure((StdNormal(), StdUniform()^2))
        g(v) = logdensityof(pushfwd(flat, Pt), v)
        v = [0.3, 0.2, 0.7]
        @test Zygote.gradient(g, v)[1] ≈ [-0.3, 0.0, 0.0]
    end
end


# Doubling that counts its evaluations:
const _doubling_calls = Ref(0)
struct CountedDoubling end
struct CountedHalving end
(::CountedDoubling)(x) = (_doubling_calls[] += 1; 2 .* x)
(::CountedHalving)(y) = y ./ 2
InverseFunctions.inverse(::CountedDoubling) = CountedHalving()
InverseFunctions.inverse(::CountedHalving) = CountedDoubling()
ChangesOfVariables.with_logabsdet_jacobian(f::CountedHalving, y) = (f(y), -length(y) * log(2))

@testset "variate sizes of pushforwards" begin
    using MeasureBase: mspace_ndims, mspace_flatsize, fixed_stream_size, restrict, logdensities, StdNormal
    using AffineMaps: Mul
    using StaticArrays: SMatrix, Size
    using Static: static

    # Origins without random variates: construction and point densities
    # need no variate size, batches take it from the type `f` returns or
    # from a declaration:
    r = restrict(x -> x > 0, StdNormal())
    ν = pushfwd(Mul(2.0), r)
    @test logdensityof(ν, 1.0) ≈ logdensityof(StdNormal(), 0.5) - log(2)
    @test logdensityof(ν, -1.0) == -Inf
    @test @inferred(mspace_ndims(typeof(ν))) == 0
    @test logdensities(ν, [1.0, -1.0, 3.0]) ≈ [logdensityof(ν, y) for y in (1.0, -1.0, 3.0)]
    r2 = restrict(x -> all(>(0), x), StdNormal()^2)
    ν2 = pushfwd(CountedDoubling(), r2)
    @test logdensityof(ν2, [1.0, 2.0]) ≈ logdensityof(StdNormal()^2, [0.5, 1.0]) - 2 * log(2)
    @test_throws ArgumentError logdensities(ν2, [1.0 2.0; 2.0 1.0])
    ν2d = pushfwd(CountedDoubling(), r2; varsize = (2,))
    @test @inferred(mspace_ndims(typeof(ν2d))) == 1 && fixed_stream_size(typeof(ν2d)) == true
    @test logdensities(ν2d, [1.0 2.0; 2.0 -1.0]) ≈ [logdensityof(ν2, [1.0, 2.0]), -Inf]
    @test_throws ArgumentError pushfwd(Mul(2.0), r; varsize = 3)

    # Sizes are learned on demand, the function is not evaluated at
    # construction or for point densities:
    _doubling_calls[] = 0
    μ3 = StdNormal()^3
    ν3 = pushfwd(CountedDoubling(), μ3)
    y = [0.2, 0.4, 0.6]
    @test logdensityof(ν3, y) ≈ logdensityof(μ3, y ./ 2) - 3 * log(2)
    @test _doubling_calls[] == 0
    @test @inferred(mspace_ndims(typeof(ν3))) == 1
    @test @inferred(mspace_flatsize(ν3)) == (3,)
    @test _doubling_calls[] == 1
    Y = [0.2 0.4; 0.4 -0.8; 0.6 1.0]
    @test logdensities(ν3, Y) ≈ [logdensityof(ν3, Y[:, j]) for j in 1:2]
    @test_throws ArgumentError logdensities(ν3, Y[1:2, :])

    # Static sizes follow from the type `f` returns for statically sized
    # variates:
    n_calls = _doubling_calls[]
    A = SMatrix{2,2}(1.0, 0.5, 0.0, 2.0)
    νs = pushfwd(Mul(A), StdNormal()^static(2))
    @test @inferred(mspace_flatsize(typeof(νs))) === Size(2)
    @test @inferred(mspace_flatsize(typeof(pushfwd(CountedDoubling(), μ3; varsize = (static(3),))))) === Size(3)
    @test @inferred(mspace_flatsize(typeof(pushfwd(CountedDoubling(), StdNormal()^static(3))))) === Size(3)
    @test _doubling_calls[] == n_calls
end
