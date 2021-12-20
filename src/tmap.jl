# A type-level map

@generated function tmap(f::F, ::Type{T}) where {F,T<:Tuple}
    f = functioninstance(F)
    Tuple{map(f, T.types)...}
end

function tmap(f::F, ::Type{A}) where {F,T,A<:AbstractArray{T}}
    f = functioninstance(F)
    p = Tuple(A.parameters)
    C = constructorof(A)
    C{(f(T), Base.tail(p)...)...}
end

function tmap(f::F, ::Type{ProductMeasure{T}}) where {F,T}
    f = functioninstance(F)
    ProductMeasure{tmap(f, T)}
end
