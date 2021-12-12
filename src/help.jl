export help_tbasemeasure_type

function help_tbasemeasure_type(μ::M) where M
    b = basemeasure(μ)
    B = typeof(b)
    println("""
    function tbasemeasure_type(::Type{$M}) 
        $B
    end

    """
    )
    return B
end
