# This file is a part of DistributionMeasures.jl, licensed under the MIT License (MIT).


MeasureBase.getdof(ν::MvNormal) = length(ν)

MeasureBase.transport_origin(ν::MvNormal) = StandardDist{Normal}(length(ν))

function MeasureBase.from_origin(ν::MvNormal, x)
    A = cholesky(ν.Σ).L
    b = ν.μ
    muladd(A, x, b)
end

function MeasureBase.to_origin(ν::MvNormal, y)
    A = cholesky(ν.Σ).L
    b = ν.μ
    A \ (y - b)
end
