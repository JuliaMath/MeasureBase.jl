# A type-level map

@inline function tmap(f::F, ::Type{T}) where {F,T<:Tuple}
    Tuple{map(f, T.types)...}
end

@inline function tmap(f::F, ::Type{A}) where {F,T,A<:AbstractArray{T}}
    p = Tuple(A.parameters)
    C = constructorof(A)
    C{(f(T), Base.tail(p)...)...}
end

@inline function tmap(f::F, ::Type{Base.Generator{G,I}}) where {F,G,I}
    Base.Generator{ComposedFunction{F,G}, I}
end

function tmap(::typeof(tbasemeasure_type), ::Type{ProductMeasure{T}}) where {T}
    if @generated
        B = tmap(tbasemeasure_type, T)
        ProductMeasure{B}
    else
        B = tmap(tbasemeasure_type, T)
        ProductMeasure{B}
    end
end
