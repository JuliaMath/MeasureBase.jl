# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase:
    weightedmeasure, superpose, productmeasure, powermeasure, pushfwd, pullbck, restrict
using MeasureBase:
    WeightedMeasure,
    SuperpositionMeasure,
    ProductMeasure,
    PowerMeasure,
    PushforwardMeasure,
    RestrictedMeasure,
    Dirac,
    StdNormal,
    StdUniform,
    StdExponential,
    PushfwdRootMeasure,
    AdaptRootMeasure
using FillArrays: Fill
using Static: static

@testset "smart constructors" begin
    @testset "powermeasure" begin
        @test powermeasure(StdNormal(), ()) === StdNormal()

        @test powermeasure(Dirac(4.2), (3,)) == Dirac(Fill(4.2, 3))

        wpw = weightedmeasure(0.3, StdNormal())^(2, 3)
        @test wpw isa WeightedMeasure
        @test wpw.logweight ≈ 6 * 0.3
        @test wpw.base == StdNormal()^(2, 3)

        # Weight pull-out must not depend on the collapsed base type:
        wd = weightedmeasure(0.3, Dirac(1.5))^3
        @test wd isa WeightedMeasure
        @test wd.logweight ≈ 3 * 0.3
        @test wd.base == Dirac(Fill(1.5, 3))

        ws = weightedmeasure(static(0.5), StdNormal())^static(4)
        @test ws.logweight ≈ 2.0
    end

    @testset "productmeasure" begin
        @test productmeasure(Fill(StdUniform(), 3)) == StdUniform()^3

        @test productmeasure(()) === Dirac(())
        @test productmeasure(NamedTuple()) === Dirac(NamedTuple())

        @test productmeasure((Dirac(1), Dirac(2))) === Dirac((1, 2))
        @test productmeasure((a = Dirac(1), b = Dirac(2))) === Dirac((a = 1, b = 2))
        @test productmeasure([Dirac(1), Dirac(2)]) == Dirac([1, 2])

        pt = productmeasure((2.0 * StdNormal(), 3.0 * StdUniform()))
        @test pt isa WeightedMeasure
        @test exp(pt.logweight) ≈ 6
        @test pt.base == ProductMeasure((StdNormal(), StdUniform()))

        pnt = productmeasure((a = 2.0 * StdNormal(), b = 3.0 * StdUniform()))
        @test pnt isa WeightedMeasure
        @test exp(pnt.logweight) ≈ 6
        @test pnt.base == ProductMeasure((a = StdNormal(), b = StdUniform()))

        pa = productmeasure([2.0 * StdNormal(), 3.0 * StdNormal()])
        @test pa isa WeightedMeasure
        @test exp(pa.logweight) ≈ 6
        @test pa.base == StdNormal()^2

        @test logdensityof(pt, (0.3, 0.5)) ≈ log(6) + logdensityof(StdNormal(), 0.3)

        @test productmeasure([StdNormal(), StdNormal()]) == StdNormal()^2
        @test productmeasure([StdNormal()^2, StdUniform()^3]) isa ProductMeasure
    end

    @testset "superpose" begin
        μ, ν = StdNormal(), StdUniform()

        @test superpose(μ) === μ
        @test superpose(μ, ν) == SuperpositionMeasure((μ, ν))

        s2 = superpose(μ, μ)
        @test s2 isa WeightedMeasure && exp(s2.logweight) ≈ 2 && s2.base === μ

        s4 = superpose(μ, μ, μ, μ)
        @test s4 isa WeightedMeasure && exp(s4.logweight) ≈ 4

        c = superpose(2.0 * μ, 3.0 * μ)
        @test c isa WeightedMeasure && exp(c.logweight) ≈ 5 && c.base === μ
        @test exp(superpose(2.0 * μ, μ).logweight) ≈ 3
        @test exp(superpose(μ, 2.0 * μ).logweight) ≈ 3
        @test superpose(2.0 * μ, 3.0 * ν) == SuperpositionMeasure((2.0 * μ, 3.0 * ν))

        ss = superpose(superpose((μ, ν)), superpose((ν, StdExponential())))
        @test ss.components === (μ, ν, ν, StdExponential())
        @test superpose(superpose((μ, ν)), StdExponential()).components ===
              (μ, ν, StdExponential())
        @test superpose(StdExponential(), superpose((μ, ν))).components ===
              (StdExponential(), μ, ν)

        # Merging must not mutate existing superpositions:
        sv = superpose(AbstractMeasure[μ, ν])
        sv2 = superpose(sv, StdExponential())
        @test length(sv.components) == 2 && length(sv2.components) == 3

        # Simplifications must be type stable, so they only happen when
        # measure equality is decidable from the measure types:
        @inferred superpose(μ, ν)
        @inferred superpose(μ, μ)
        @inferred superpose(2.0 * μ, 3.0 * μ)
        @inferred superpose(2.0 * μ, μ)
        @inferred superpose(μ, μ, μ, μ)
        @inferred superpose(Dirac(1), Dirac(1))
        @test superpose(Dirac(1), Dirac(1)) isa SuperpositionMeasure
        @inferred superpose(2.0 * Dirac(1), 3.0 * Dirac(1))
        @test superpose(2.0 * Dirac(1), 3.0 * Dirac(1)) isa SuperpositionMeasure

        @test superpose(Fill(μ, 4)) == weightedmeasure(log(4), μ)
        @test superpose([μ, μ, μ]) == weightedmeasure(log(3), μ)

        @test logdensityof(c, 0.3) ≈ log(5) + logdensityof(μ, 0.3)
    end

    @testset "pushfwd" begin
        @test pushfwd(identity, StdNormal()) === StdNormal()
        @test pushfwd(identity, StdNormal(), PushfwdRootMeasure()) === StdNormal()

        @test pushfwd(sqrt, Dirac(4.0)) === Dirac(2.0)
        @test pushfwd(sqrt, Dirac(4.0), PushfwdRootMeasure()) === Dirac(2.0)

        pw = pushfwd(sqrt, 3.0 * StdExponential())
        @test pw isa WeightedMeasure && exp(pw.logweight) ≈ 3
        @test pw.base isa PushforwardMeasure

        pp = pushfwd(exp, pushfwd(sqrt, StdExponential()))
        @test pp isa PushforwardMeasure && pp.origin === StdExponential()

        @test pullbck(log, Dirac(4.0)) === Dirac(exp(4.0))
    end

    @testset "restrict" begin
        r = restrict(x -> x > 0, StdNormal())
        @test r isa RestrictedMeasure && r.base === StdNormal()

        r2 = restrict(x -> x < 1, r)
        @test r2 isa RestrictedMeasure && r2.base === StdNormal()
        @test r2.predicate(0.5) && !r2.predicate(-1.0) && !r2.predicate(2.0)

        @test restrict(x -> x > 0)(StdNormal()) isa RestrictedMeasure
    end
end
