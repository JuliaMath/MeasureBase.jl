smf_measures = [
    StdNormal()
    StdLogistic()
    StdUniform()
    Half(StdNormal())
    Half(StdLogistic())
]

@testset "smf" begin
    for μ in smf_measures
        test_smf(μ)
    end
end
