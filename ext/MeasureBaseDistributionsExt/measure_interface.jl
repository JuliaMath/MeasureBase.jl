# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Bare distributions behave like their measures (some constructors keep
# them, e.g. products of distributions). Their densities are Distributions'
# ones, relative to the base measures of their measures:
@inline MeasureBase.logdensity_def(d::Distribution, x) = DensityInterface.logdensityof(d, x)
@inline MeasureBase.unsafe_logdensityof(d::Distribution, x) = DensityInterface.logdensityof(d, x)
@inline MeasureBase.insupportof(d::Distribution, x) = Distributions.insupport(d, x) & _finite_variate(d, x)
@inline MeasureBase.basemeasure(d::Distribution) = basemeasure(asmeasure(d))
@inline MeasureBase.rootmeasure(d::Distribution) = rootmeasure(asmeasure(d))
@inline MeasureBase.testvalue(d::Distribution) = testvalue(asmeasure(d))
@inline MeasureBase.testvalue(::Type{T}, d::Distribution) where {T} = testvalue(T, asmeasure(d))

MeasureBase.mintegrate(f, base::Distribution) =
    MeasureBase.mintegrate(f, convert(AbstractMeasure, base))


# Parameters as returned by `Distributions.params`, named after the
# leading fields of the distribution where these hold them (by type),
# after the constructor arguments for some common distributions whose
# fields differ, positionally (`p1`, `p2`, ...) otherwise:
function MeasureBase.paramsof(d::Distribution)
    p = Distributions.params(d)
    NamedTuple{_paramnames(typeof(d), typeof(p))}(p)
end
MeasureBase.paramnames(d::Distribution) = keys(MeasureBase.paramsof(d))

@generated function _paramnames(::Type{D}, ::Type{P}) where {D<:Distribution,P<:Tuple}
    n = fieldcount(P)
    by_fields = fieldcount(D) >= n && all(i -> fieldtype(D, i) == fieldtype(P, i), 1:n)
    by_fields ? fieldnames(D)[1:n] : ntuple(i -> Symbol(:p, i), n)
end
_paramnames(::Type{<:Distributions.Categorical}, ::Type{<:Tuple{Any}}) = (:p,)
_paramnames(::Type{<:Distributions.AbstractMvTDist}, ::Type{<:NTuple{3,Any}}) = (:df, :μ, :Σ)
