export help_basemeasure_type

function help_basemeasure_type(μ::M) where M
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
