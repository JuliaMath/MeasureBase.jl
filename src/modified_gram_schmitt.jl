using LinearAlgebra

# Modify column jᵥ to be orthogonal to column jᵤ
function ortho!(A::AbstractMatrix, jᵥ, jᵤ; startat=1)
    # v ← v - projᵤ(v)
    nrows = size(A,1)
    u_v = 0.0
    u_u = 0.0
    @inbounds for i in startat:nrows
        u_v += A[i,jᵤ] * A[i,jᵥ]
        u_u += A[i,jᵤ] ^ 2
    end

    k = u_v / u_u

    @inbounds for i in startat:nrows
        A[i, jᵥ] = A[i, jᵥ] - k * A[i, jᵤ]
    end

    return A
end

function ortho!(A::AbstractMatrix, triangular=false)
    for j2 in size(A,2):-1:1
        for j1 in 1:(j2-1)
            ortho!(A, j1, j2)
        end
    end

    for col in eachcol(A)
        normalize!(col)
    end

    return A
end

function makelower!(A)
    for jᵥ in 2:size(A,2)
        @inbounds for jᵤ in 1:(jᵥ-1)
            k = A[jᵤ, jᵥ] / A[jᵤ, jᵤ]

            A[jᵤ,jᵥ] = 0.0
            
            @inbounds for i in (jᵤ+1):size(A,1)
                A[i,jᵥ] = A[i, jᵥ] - k * A[i, jᵤ]
            end
        end
    end

    for j in 1:size(A,2)
        @inbounds if A[j,j] < 0.0
            view(A, :, j) .*= -1
        end
    end

    return A

end




function WtoV(W)
    # make columns orthogonal
    for j2 in size(W,2):-1:1
        for j1 in 1:(j2-1)
            ortho!(W, j1, j2; startat=j2)
        end
    end

    for col in eachcol(W)
        normalize!(col)
    end

    return W
end

W = [I; randn(5,5) ./ randn(5,5)]

using BenchmarkTools
@btime WtoV(W)



W' * W