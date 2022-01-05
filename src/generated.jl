@generated function _logdensityof(μ, x, ℓ, ::StaticInt{n}) where {n}
    nsteps = max(n, 0)
    quote
        $(Expr(:meta, :inline))
        Base.Cartesian.@nexprs $nsteps i -> begin
            Δℓ = logdensity_def(μ, x)
            # @show μ
            # @show Δℓ
            # println()
            μ = basemeasure(μ, x)
            ℓ += Δℓ
        end
        return (ℓ,μ)
    end
end

@generated function _rootmeasure(μ, ::StaticInt{n}) where {n}
    q = quote end
    foreach(1:n) do _
        push!(q.args, :(μ = basemeasure(μ)))
    end
    return q
end

