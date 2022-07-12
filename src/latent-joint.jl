latentof(m) = m
manifestof(m) = m

function jointof(m)
    fwd(x) = x => x
    
    function back(p::Pair)
        x,y = p
        @assert x === y
        return x
    end

    PushforwardMeasure(fwd, back, m, NoVolCorr())
end
