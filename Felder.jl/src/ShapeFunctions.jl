export AbstractShapeFunctions
export AbstractScalarShapeFunctions
export ConstantShapeFunction
export getorder
export ndofs
export ndofs_vertex, ndofs_edge, ndofs_face, ndofs_cell
export ndofs_facet, ndofs_element
export sftype
export shapefunc_N
export shapefunc_∇N
export shapefunc_∇∇N
export eachshapefunc_N
export eachshapefunc_∇N
export eachshapefunc_∇∇N
export eachshapefunc_N_∇N
export refcoordinates_grid

abstract type AbstractShapeFunctions{Order} end

getorder(::Union{T, Type{T}}) where {T <: AbstractShapeFunctions{N}} where {N} = N

abstract type AbstractScalarShapeFunctions{Order} <: AbstractShapeFunctions{Order} end
abstract type AbstractVectorShapeFunctions{Order} <: AbstractShapeFunctions{Order} end

function each_dof_refcoordinates(::Type{S}, sf::AbstractShapeFunctions) where {S<:AbstractShape}
    (dof_refcoordinates(S, sf, i) for i in 1:ndofs(S, sf))
end

"""
"""
function eachshapefunc_N(::Type{S}, sf::AbstractShapeFunctions, ξ::AbstractVector) where {S<:AbstractShape}
    (shapefunc_N(S, sf, i, ξ) for i in 1:ndofs(S, sf))
end

"""
"""
function eachshapefunc_∇N(::Type{S}, sf::AbstractScalarShapeFunctions, ξ::AbstractVector) where {S<:AbstractShape}
    (shapefunc_∇N(S, sf, i, ξ) for i in 1:ndofs(S, sf))
end

"""
"""
function eachshapefunc_N_∇N(::Type{S}, sf::AbstractScalarShapeFunctions, ξ::AbstractVector) where {S<:AbstractShape}
    ((shapefunc_N(S, sf, i, ξ), shapefunc_∇N(S, sf, i, ξ)) for i in 1:ndofs(S, sf))
end

"""
"""
function eachshapefunc_∇∇N(::Type{S}, sf::AbstractScalarShapeFunctions, ξ::AbstractVector) where {S<:AbstractShape}
    (shapefunc_∇∇N(S, sf, i, ξ) for i in 1:ndofs(S, sf))
end

##############################
# Constant Shape Functions
##############################

struct ConstantShapeFunction <: AbstractScalarShapeFunctions{0} end

ndofs(::Type{<:AbstractShape}, ::ConstantShapeFunction) = 1

dof_refcoordinates(::Type{S}, ::ConstantShapeFunction, ::Integer) where {S <: AbstractShape} = center_refcoordinates(S)

function shapefunc_N(::Type{S}, ::ConstantShapeFunction, i::Integer, λ::AbstractVector) where {S <: AbstractShape}
    return 1.0
end

function shapefunc_∇N(::Type{S}, ::ConstantShapeFunction, i::Integer, λ::AbstractVector) where {S <: AbstractShape}
    return zero(SVector{ndims(S), Float64})
end

edgedofs(::Type{<:AbstractShape}, ::ConstantShapeFunction, i::Integer) = SA[0]
facetdofs(::Type{<:AbstractShape}, ::ConstantShapeFunction, i::Integer) = SA[0]
