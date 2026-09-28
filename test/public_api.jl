using Test
using MeasureBase

@testset "public API" begin
    @static if VERSION >= v"1.11"
        public_names = filter(n -> Base.ispublic(MeasureBase, n), names(MeasureBase; all = true))
        @test all(n -> isdefined(MeasureBase, n), public_names)
        # Extension points and internals used by downstream packages:
        for n in (
            :logdensityof_impl, :batched_logdensityof_impl, :rand_impl, :batched_rand_impl,
            :transport_to_std, :transport_from_std, :preferred_stdmeasure, :mspace_ndims,
            :mspace_flatsize, :checked_arg, :NoDOF, :UnknownMass, :DensityMeasure,
            :TransportFunction, :GenContext,
        )
            @test Base.ispublic(MeasureBase, n) && !Base.isexported(MeasureBase, n)
        end
    end
end
