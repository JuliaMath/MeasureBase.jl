diagonal_pair(x) = x => x

latentof(m::AbstractMeasure) = m
manifestof(m::AbstractMeasure) = m
jointof(m::AbstractMeasure) = pushfwd(diagonal_pair, m)