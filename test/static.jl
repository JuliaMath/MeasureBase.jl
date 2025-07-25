using Test

import MeasureBase
using MeasureBase:
    StaticUnitRange,
    StaticOneTo,
    IntegerLike,
    SizeLike,
    StaticSizeLike,
    AxesLike,
    StaticAxesLike,
    OneToLike,
    StaticOneToLike,
    StaticUnitRangeLike,
    one_to,
    asnonstatic,
    fill_with,
    staticarray_type,
    maybestatic_reshape,
    maybestatic_length,
    maybestatic_size,
    maybestatic_axes,
    axes2size,
    size2axes,
    size2length,
    asaxes,
    maybestatic_eachindex,
    maybestatic_first,
    maybestatic_last,
    NoTypeSize,
    size_from_type,
    element_size,
    canonical_indices,
    canonical_size,
    canonical_axes

import Static
using Static: static
import StaticArrays
import FillArrays

@testset "static" begin
    v = 4.2
    T = typeof(v)

    tpl = (7, 42, 5)
    nt = (a = 7, b = 42, c = 5)

    i = 7
    si = static(7)

    @test i isa IntegerLike
    @test si isa IntegerLike

    sz = (2, 4, 3)
    sasz = StaticArrays.Size(2, 4, 3)
    sisz = (static(2), static(4), static(3))

    len = prod(sz)
    slen = static(len)

    @test sz isa SizeLike
    @test sasz isa SizeLike
    @test sisz isa SizeLike

    @test !(sz isa StaticSizeLike)
    @test sasz isa StaticSizeLike
    @test sisz isa StaticSizeLike

    axs = (Base.OneTo(2), 2:5, Base.OneTo(3))
    axs1 = (Base.OneTo(2), Base.OneTo(4), Base.OneTo(3))
    saaxs = (StaticOneTo(2), StaticUnitRange(2, 5), StaticOneTo(3))
    saaxs1 = (StaticOneTo(2), StaticOneTo(4), StaticOneTo(3))
    siaxs = (Static.SOneTo(2), static(2):static(5), static(1):static(3))
    siaxs1 = (Static.SOneTo(2), static(1):static(4), static(1):static(3))

    @test axs isa AxesLike
    @test axs1 isa AxesLike
    @test saaxs isa AxesLike
    @test saaxs1 isa AxesLike
    @test siaxs isa AxesLike
    @test siaxs1 isa AxesLike

    @test !(axs isa StaticAxesLike)
    @test saaxs isa StaticAxesLike
    @test saaxs1 isa StaticAxesLike
    @test siaxs isa StaticAxesLike
    @test siaxs1 isa StaticAxesLike

    @test axs[1] isa OneToLike
    @test !(axs[2] isa OneToLike)
    @test axs[3] isa OneToLike

    @test saaxs[1] isa OneToLike
    @test !(saaxs[2] isa OneToLike)
    @test saaxs1[2] isa OneToLike
    @test saaxs[3] isa OneToLike

    @test siaxs[1] isa OneToLike
    @test !(siaxs[2] isa OneToLike)
    @test siaxs1[2] isa OneToLike
    @test siaxs[3] isa OneToLike

    @test !(axs[1] isa StaticOneToLike)
    @test !(axs[2] isa StaticOneToLike)
    @test !(axs[3] isa StaticOneToLike)

    @test saaxs[1] isa StaticOneToLike
    @test !(saaxs[2] isa StaticOneToLike)
    @test saaxs1[2] isa StaticOneToLike
    @test saaxs[3] isa StaticOneToLike

    @test siaxs[1] isa StaticOneToLike
    @test !(siaxs[2] isa StaticOneToLike)
    @test siaxs1[2] isa StaticOneToLike
    @test siaxs[3] isa StaticOneToLike

    @test !(axs[1] isa StaticUnitRangeLike)
    @test !(axs[2] isa StaticUnitRangeLike)
    @test !(axs[3] isa StaticUnitRangeLike)

    @test saaxs[1] isa StaticUnitRangeLike
    @test saaxs[2] isa StaticUnitRangeLike
    @test saaxs1[2] isa StaticUnitRangeLike
    @test saaxs[3] isa StaticUnitRangeLike

    @test siaxs[1] isa StaticUnitRangeLike
    @test siaxs[2] isa StaticUnitRangeLike
    @test siaxs1[2] isa StaticUnitRangeLike
    @test siaxs[3] isa StaticUnitRangeLike

    @test @inferred(one_to(i)) == Base.OneTo(i)
    @test @inferred(one_to(si)) == StaticOneTo(i)

    @test @inferred(asnonstatic(i)) === i
    @test @inferred(asnonstatic(si)) === i
    @test @inferred(asnonstatic(sz)) === sz
    @test @inferred(asnonstatic(sasz)) === sz
    @test @inferred(asnonstatic(sisz)) === sz
    @test @inferred(asnonstatic(axs)) === axs
    @test @inferred(asnonstatic(saaxs)) === axs
    @test @inferred(asnonstatic(saaxs1)) === (Base.OneTo(2), Base.OneTo(4), Base.OneTo(3))
    @test @inferred(asnonstatic(siaxs)) === axs
    @test @inferred(asnonstatic(siaxs1)) === (Base.OneTo(2), Base.OneTo(4), Base.OneTo(3))

    @test @inferred(fill_with(v, i)) === FillArrays.Fill(v, i)
    @test @inferred(fill_with(v, si)) === StaticArrays.SVector(fill(v, i)...)
    @test @inferred(fill_with(v, ())) === FillArrays.Fill(v)

    @test @inferred(fill_with(v, sz)) === FillArrays.Fill(v, sz)
    @test @inferred(fill_with(v, sasz)) === StaticArrays.SArray{Tuple{sz...},T}(fill(v, sz))
    @test @inferred(fill_with(v, sisz)) === StaticArrays.SArray{Tuple{sz...},T}(fill(v, sz))

    @test @inferred(fill_with(v, axs)) === FillArrays.Fill(v, axs)
    @test @inferred(fill_with(v, saaxs)) === FillArrays.Fill(v, axs)
    @test @inferred(fill_with(v, saaxs1)) ===
          StaticArrays.SArray{Tuple{sz...},T}(fill(v, sz))
    @test @inferred(fill_with(v, siaxs)) === FillArrays.Fill(v, axs)
    @test @inferred(fill_with(v, siaxs1)) ===
          StaticArrays.SArray{Tuple{sz...},T}(fill(v, sz))

    @test @inferred(staticarray_type(T, sasz)) <: StaticArrays.SArray{Tuple{2,4,3},T}

    A = rand(T, len)
    FA = FillArrays.Fill(v, len)
    SA = StaticArrays.SVector(A...)

    # Array with CartesianIndices
    ciA = view(rand(5, 6, 6), 3:4, 2:5, 3:5)
    ciidxs = eachindex(ciA)

    rshpA = reshape(A, sz)
    rshpFA = FillArrays.Fill(v, sz)
    rshpSA = StaticArrays.SArray{Tuple{sz...},T}(A)

    @test @inferred(maybestatic_reshape(A, sz)) == rshpA
    @test typeof(maybestatic_reshape(A, sz)) == typeof(rshpA)
    @test @inferred(maybestatic_reshape(A, sasz)) == rshpA
    @test maybestatic_reshape(A, sasz) isa StaticArrays.SArray
    @test @inferred(maybestatic_reshape(A, sisz)) == rshpA
    @test maybestatic_reshape(A, sisz) isa StaticArrays.SArray

    @test @inferred(maybestatic_reshape(FA, sz)) == rshpFA
    @test typeof(maybestatic_reshape(FA, sz)) == typeof(rshpFA)
    @test @inferred(maybestatic_reshape(FA, sasz)) == rshpFA
    @test maybestatic_reshape(FA, sasz) isa StaticArrays.SArray
    @test @inferred(maybestatic_reshape(FA, sisz)) == rshpFA
    @test maybestatic_reshape(FA, sisz) isa StaticArrays.SArray

    @test @inferred(maybestatic_reshape(SA, sz)) == rshpA
    @test maybestatic_reshape(SA, sz) isa Base.ReshapedArray{T,3,<:StaticArrays.SVector}
    @test @inferred(maybestatic_reshape(SA, sasz)) === rshpSA
    @test @inferred(maybestatic_reshape(SA, sisz)) === rshpSA

    @test @inferred(maybestatic_length(5)) === static(1)
    @test @inferred(maybestatic_length(())) === static(0)
    @test @inferred(maybestatic_length((sz))) === static(3)
    @test @inferred(maybestatic_length((a = 2, b = 4, c = 3))) === static(3)
    @test @inferred(maybestatic_length(Base.OneTo(4))) === 4
    @test @inferred(maybestatic_length(StaticArrays.SOneTo(4))) === static(4)
    @test @inferred(maybestatic_length(Static.SOneTo(4))) === static(4)
    @test @inferred(maybestatic_length(static(2):static(5))) === static(4)
    @test @inferred(maybestatic_length(rshpA)) === length(rshpA)
    @test @inferred(maybestatic_length(rshpFA)) === length(rshpA)
    @test @inferred(maybestatic_length(rshpSA)) === static(length(rshpA))

    @test @inferred(maybestatic_size(5)) === ()
    @test @inferred(maybestatic_size(())) === StaticArrays.Size(0)
    @test @inferred(maybestatic_size((sz))) === StaticArrays.Size(3)
    @test @inferred(maybestatic_size((sasz))) === StaticArrays.Size(3)
    @test @inferred(maybestatic_size((sisz))) === StaticArrays.Size(3)
    @test @inferred(maybestatic_size((a = 2, b = 4, c = 3))) === StaticArrays.Size(3)
    @test @inferred(maybestatic_size(Base.OneTo(4))) === (4,)
    @test @inferred(maybestatic_size(StaticArrays.SOneTo(4))) === StaticArrays.Size(4)
    @test @inferred(maybestatic_size(StaticUnitRange(2, 5))) === StaticArrays.Size(4)
    @test @inferred(maybestatic_size(Static.SOneTo(4))) === StaticArrays.Size(4)
    @test @inferred(maybestatic_size(static(2):static(5))) === StaticArrays.Size(4)
    @test @inferred(maybestatic_size(rshpA)) === size(rshpA)
    @test @inferred(maybestatic_size(rshpFA)) === size(rshpA)
    @test @inferred(maybestatic_size(rshpSA)) === StaticArrays.Size(size(rshpA)...)

    @test @inferred(maybestatic_axes(5)) === ()
    @test @inferred(maybestatic_axes(())) === (StaticOneTo(0),)
    @test @inferred(maybestatic_axes((sz))) === (StaticOneTo(3),)
    @test @inferred(maybestatic_axes((sasz))) === (StaticOneTo(3),)
    @test @inferred(maybestatic_axes((sisz))) === (StaticOneTo(3),)
    @test @inferred(maybestatic_axes((a = 2, b = 4, c = 3))) === (StaticOneTo(3),)
    @test @inferred(maybestatic_axes(Base.OneTo(4))) === (Base.OneTo(4),)
    @test @inferred(maybestatic_axes(StaticArrays.SOneTo(4))) === (StaticOneTo(4),)
    @test @inferred(maybestatic_axes(Static.SOneTo(4))) === (StaticOneTo(4),)
    @test @inferred(maybestatic_axes(static(2):static(5))) === (StaticOneTo(4),)
    @test @inferred(maybestatic_axes(rshpA)) === axes(rshpA)
    @test @inferred(maybestatic_axes(rshpFA)) === axes(rshpA)
    @test @inferred(maybestatic_axes(rshpSA)) === saaxs1

    @test @inferred(axes2size(())) === ()
    @test @inferred(axes2size(axs)) === sz
    @test @inferred(axes2size(saaxs)) === sasz
    @test @inferred(axes2size(saaxs1)) === sasz
    @test @inferred(axes2size(siaxs)) === sasz
    @test @inferred(axes2size(siaxs1)) === sasz

    @test @inferred(size2axes(())) === ()
    @test @inferred(size2axes(sz)) === axs1
    @test @inferred(size2axes(sasz)) === saaxs1
    @test @inferred(size2axes(sisz)) === saaxs1

    @test @inferred(size2length(())) === static(1)
    @test @inferred(size2length(sz)) === len
    @test @inferred(size2length(sasz)) === slen
    @test @inferred(size2length(sisz)) === slen

    @test @inferred(asaxes(())) === ()
    @test @inferred(asaxes(len)) === (Base.OneTo(len),)
    @test @inferred(asaxes(slen)) === (StaticOneTo(len),)
    @test @inferred(asaxes(sz)) === axs1
    @test @inferred(asaxes(sasz)) === saaxs1
    @test @inferred(asaxes(sisz)) === saaxs1
    @test @inferred(asaxes(axs)) === axs
    @test @inferred(asaxes(axs1)) === axs1
    @test @inferred(asaxes(saaxs)) === saaxs
    @test @inferred(asaxes(saaxs1)) === saaxs1
    @test @inferred(asaxes(siaxs)) === siaxs
    @test @inferred(asaxes(siaxs1)) === siaxs1

    @test @inferred(maybestatic_eachindex(())) === StaticOneTo(0)
    @test @inferred(maybestatic_eachindex(tpl)) === StaticOneTo(3)
    @test @inferred(maybestatic_eachindex(nt)) === StaticOneTo(3)
    @test @inferred(maybestatic_eachindex(axs[1])) === Base.OneTo(length(axs[1]))
    @test @inferred(maybestatic_eachindex(axs[2])) === Base.OneTo(length(axs[2]))
    @test @inferred(maybestatic_eachindex(saaxs[1])) === StaticOneTo(length(axs[1]))
    @test @inferred(maybestatic_eachindex(saaxs[2])) === StaticOneTo(length(axs[2]))
    @test @inferred(maybestatic_eachindex(siaxs[1])) === StaticOneTo(length(axs[1]))
    @test @inferred(maybestatic_eachindex(siaxs[2])) === StaticOneTo(length(axs[2]))
    @test @inferred(maybestatic_eachindex(A)) === Base.OneTo(24)
    @test @inferred(maybestatic_eachindex(ciA)) === eachindex(ciA)
    @test @inferred(maybestatic_eachindex(FA)) === Base.OneTo(24)
    @test @inferred(maybestatic_eachindex(SA)) === StaticOneTo(24)

    @test_throws BoundsError maybestatic_first(())
    @test @inferred(maybestatic_first(tpl)) === first(tpl)
    @test @inferred(maybestatic_first(nt)) === first(nt)
    @test @inferred(maybestatic_first(sz)) === first(sz)
    @test @inferred(maybestatic_first(sasz)) === static(first(sz))
    @test @inferred(maybestatic_first(sisz)) === static(first(sz))
    @test @inferred(maybestatic_first(axs[1])) === first(axs[1])
    @test @inferred(maybestatic_first(axs[2])) === first(axs[2])
    @test @inferred(maybestatic_first(saaxs[1])) === static(first(axs[1]))
    @test @inferred(maybestatic_first(saaxs[2])) === static(first(axs[2]))
    @test @inferred(maybestatic_first(siaxs[1])) === static(first(axs[1]))
    @test @inferred(maybestatic_first(siaxs[2])) === static(first(axs[2]))
    @test @inferred(maybestatic_first(A)) === first(A)
    @test @inferred(maybestatic_first(ciA)) === first(ciA)
    @test @inferred(maybestatic_first(FA)) === first(FA)
    @test @inferred(maybestatic_first(SA)) === first(SA)

    @test_throws BoundsError maybestatic_last(())
    @test @inferred(maybestatic_last(tpl)) === last(tpl)
    @test @inferred(maybestatic_last(nt)) === last(nt)
    @test @inferred(maybestatic_last(sz)) === last(sz)
    @test @inferred(maybestatic_last(sasz)) === static(last(sz))
    @test @inferred(maybestatic_last(sisz)) === static(last(sz))
    @test @inferred(maybestatic_last(axs[1])) === last(axs[1])
    @test @inferred(maybestatic_last(axs[2])) === last(axs[2])
    @test @inferred(maybestatic_last(saaxs[1])) === static(last(axs[1]))
    @test @inferred(maybestatic_last(saaxs[2])) === static(last(axs[2]))
    @test @inferred(maybestatic_last(siaxs[1])) === static(last(axs[1]))
    @test @inferred(maybestatic_last(siaxs[2])) === static(last(axs[2]))
    @test @inferred(maybestatic_last(A)) === last(A)
    @test @inferred(maybestatic_last(ciA)) === last(ciA)
    @test @inferred(maybestatic_last(FA)) === last(FA)
    @test @inferred(maybestatic_last(SA)) === last(SA)

    @test @inferred(canonical_indices(axs[1])) === axs[1]
    @test @inferred(canonical_indices(axs[2])) === axs[2]
    @test @inferred(canonical_indices(saaxs[1])) === saaxs[1]
    @test @inferred(canonical_indices(saaxs[2])) === saaxs[2]
    @test @inferred(canonical_indices(siaxs[1])) === saaxs[1]
    @test @inferred(canonical_indices(siaxs[2])) === saaxs[2]
    @test @inferred(canonical_indices(ciidxs)) === ciidxs

    @test @inferred(canonical_size(sz)) === sz
    @test @inferred(canonical_size(sasz)) === sasz
    @test @inferred(canonical_size(sisz)) === sasz

    @test @inferred(canonical_axes(axs)) === axs
    @test @inferred(canonical_axes(axs1)) === axs1
    @test @inferred(canonical_axes(saaxs)) === saaxs
    @test @inferred(canonical_axes(saaxs1)) === saaxs1
    @test @inferred(canonical_axes(siaxs)) === saaxs
    @test @inferred(canonical_axes(siaxs1)) === saaxs1

    @test @inferred(size_from_type(typeof(i))) === maybestatic_size(i)
    @test @inferred(size_from_type(typeof(()))) === maybestatic_size(())
    @test @inferred(size_from_type(typeof(tpl))) === maybestatic_size(tpl)
    @test @inferred(size_from_type(typeof(nt))) === maybestatic_size(nt)
    @test @inferred(size_from_type(typeof(sz))) === maybestatic_size(sz)
    @test @inferred(size_from_type(typeof(sasz))) === maybestatic_size(sasz)
    @test @inferred(size_from_type(typeof(axs))) === maybestatic_size(axs)
    @test @inferred(size_from_type(typeof(axs1))) === maybestatic_size(axs1)
    @test @inferred(size_from_type(typeof(saaxs))) === maybestatic_size(saaxs)
    @test @inferred(size_from_type(typeof(saaxs1))) === maybestatic_size(saaxs1)
    @test @inferred(size_from_type(typeof(siaxs))) === maybestatic_size(siaxs)
    @test @inferred(size_from_type(eltype(A))) === maybestatic_size(A[1])
    @test @inferred(size_from_type(typeof(A))) === NoTypeSize{Vector{T}}()

    @test @inferred(element_size(tpl)) === maybestatic_size(tpl[1])
    @test @inferred(element_size(nt)) === maybestatic_size(nt[1])
    @test @inferred(element_size(sz)) === maybestatic_size(sz[1])
    @test @inferred(element_size(sasz)) === maybestatic_size(sasz[1])
    @test @inferred(element_size(axs)) === NoTypeSize #!!!!!!!!
    @test @inferred(element_size(saaxs)) === NoTypeSize #!!!!!!!!
    @test @inferred(element_size(siaxs)) === NoTypeSize #!!!!!!!!

    @test @inferred(element_size((sz, sz))) === maybestatic_size(sz)
    @test @inferred(element_size((sasz, sasz))) === maybestatic_size(sasz)
    @test @inferred(element_size((sisz, sisz))) === maybestatic_size(sisz)
    @test @inferred(element_size((sz, sasz))) === maybestatic_size(sasz)
    @test @inferred(element_size((a = sz, b = sz))) === maybestatic_size(sz)
    @test @inferred(element_size((a = sasz, b = sasz))) === maybestatic_size(sasz)
    @test @inferred(element_size((a = sisz, b = sisz))) === maybestatic_size(sisz)
    @test @inferred(element_size((a = sz, b = sasz))) === maybestatic_size(sasz)
    @test @inferred(element_size([sz, sz])) === maybestatic_size(sz)
    @test @inferred(element_size([sasz, sasz])) === maybestatic_size(sasz)
    @test @inferred(element_size([sisz, sisz])) === maybestatic_size(sisz)
    @test @inferred(element_size([sz, sasz])) === maybestatic_size(sasz)

    @test @inferred(element_size((axs, axs))) === maybestatic_size(axs)
    @test @inferred(element_size((saaxs, saaxs))) === maybestatic_size(saaxs)
    @test @inferred(element_size((siaxs, siaxs))) === maybestatic_size(siaxs)
    @test @inferred(element_size((axs, saaxs))) === maybestatic_size(saaxs)
    @test @inferred(element_size((a = axs, b = axs))) === maybestatic_size(axs)
    @test @inferred(element_size((a = saaxs, b = saaxs))) === maybestatic_size(saaxs)
    @test @inferred(element_size((a = siaxs, b = siaxs))) === maybestatic_size(siaxs)
    @test @inferred(element_size((a = axs, b = saaxs))) === maybestatic_size(saaxs)
    @test @inferred(element_size([axs, axs])) === maybestatic_size(axs)
    @test @inferred(element_size([saaxs, saaxs])) === maybestatic_size(saaxs)
    @test @inferred(element_size([siaxs, siaxs])) === maybestatic_size(siaxs)
    @test @inferred(element_size([axs, saaxs])) === maybestatic_size(saaxs)

    @test @inferred(element_size(A)) === maybestatic_size(A[1])
    @test @inferred(element_size(fill(saaxs, 5))) === maybestatic_size(saaxs)
    @test @inferred(element_size(fill([1,2,3,4], 5))) === maybestatic_size([1,2,3,4])
    @test @inferred(element_size(Fill(saaxs, 5))) === maybestatic_size(saaxs)
    @test @inferred(element_size(Fill([1,2,3,4], 5))) === maybestatic_size([1,2,3,4])
    @test @inferred(element_size([1, 2], [3, 4, 5])) === NoTypeSize{Vector{Int}}()
end
