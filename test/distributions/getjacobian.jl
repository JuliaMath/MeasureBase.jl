# This file is a part of ChangesOfVariables.jl, licensed under the MIT License (MIT).

import ForwardDiff

torv_and_back(V::AbstractVector{<:Real}) = V, identity
torv_and_back(x::Real) = [x], V -> V[1]
torv_and_back(x::Complex) = [real(x), imag(x)], V -> Complex(V[1], V[2])
torv_and_back(x::NTuple{N}) where N = [x...], V -> ntuple(i -> V[i], Val(N))

function torv_and_back(x::Ref)
    xval = x[]
    V, to_xval = torv_and_back(xval)
    back_to_ref(V) = Ref(to_xval(V))
    return (V, back_to_ref)
end

torv_and_back(A::AbstractArray{<:Real}) = vec(A), V -> reshape(V, size(A))

function torv_and_back(A::AbstractArray{Complex{T}, N}) where {T<:Real, N}
    RA = cat(real.(A), imag.(A), dims = N+1) 
    V, to_array = torv_and_back(RA)
    function back_to_complex(V)
        RA = to_array(V)
        Complex.(view(RA, map(_ -> :, size(A))..., 1), view(RA, map(_ -> :, size(A))..., 2))
    end
    return (V, back_to_complex)
end


function getjacobian(f, x)
    V, to_x = torv_and_back(x)
    vf(V) = torv_and_back(f(to_x(V)))[1]
    ForwardDiff.jacobian(vf, V)
end
