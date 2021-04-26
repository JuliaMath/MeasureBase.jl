struct Integral{F,M}
    f::F
    μ::M
end

∫(μ) = ∫(identity,μ)
∫(f, μ) = Integral(f, μ)
