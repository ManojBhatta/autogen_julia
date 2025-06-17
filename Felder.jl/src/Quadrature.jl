export AbstractQuadrature
export NoQuadrature
export AbstractGaussQuadrature
export GaussQuadrature

export getorder
export quadweight
export quadpoint
export nquadpoints
export eachquadweight
export eachquadpoint
export eachquadweightpoint
export defaultquadrature

abstract type AbstractQuadrature{Order} end

struct NoQuadrature <: AbstractQuadrature{0} end

getorder(::Union{T, Type{T}}) where {T <: AbstractQuadrature{N}} where {N} = N

"""
"""
function eachquadweight(::Type{S}, quad::AbstractQuadrature) where {S <: AbstractShape}
    (quadweight(S, quad, i) for i in 1:nquadpoints(S, quad))
end

"""
"""
function eachquadpoint(::Type{S}, quad::AbstractQuadrature) where {S <: AbstractShape}
    (quadpoint(S, quad, i) for i in 1:nquadpoints(S, quad))
end

"""
"""
function eachquadweightpoint(::Type{S}, quad::AbstractQuadrature) where {S <: AbstractShape}
    ((quadweight(S, quad, i), quadpoint(S, quad, i)) for i in 1:nquadpoints(S, quad))
end

##############################
# Gauss Quadrature
##############################

abstract type AbstractGaussQuadrature{Order} <: AbstractQuadrature{Order} end

struct GaussQuadrature{Order} <: AbstractGaussQuadrature{Order} end

GaussQuadrature(order) = GaussQuadrature{order}()

##############################
# Vertex 0D
##############################
#=
    ○
    1
=#

nquadpoints(::Type{<:AbstractVertex}, ::GaussQuadrature) = 1

function quadweight(::Type{<:AbstractVertex}, ::GaussQuadrature, i::Integer)::Float64
    i == 1 && return 1.0
    throw(ArgumentError("no quadrature point $i for GaussQuadrature in AbstractVertex"))
end

function quadpoint(::Type{<:AbstractVertex}, ::GaussQuadrature, i::Integer)::SVector{0, Float64}
    i == 1 && return ()
    throw(ArgumentError("no quadrature point $i for GaussQuadrature in AbstractVertex"))
end

include("Gauss1D.jl")
include("Gauss2D.jl")
include("Gauss3D.jl")

"""
    defaultquadrature(shapetype)
    defaultquadrature(shapefunctiontype)
    defaultquadrature(shapefunctiontype, shapetype)

Return a default quadrature for a given shape or shape function type.

TODO: Do proper investigation The convergence with the rules below is not clear yet!
For Tri and Tet elements, a reasonable order of the quadrature rule should be 2*order of the shape function.
For Quad and Hex elements, a reasonable order of the quadrature rule should be  ....????

TODO: How to handle geometric shapefunction and field interpolation shape function types?

# Related:
https://scicomp.stackexchange.com/questions/561/are-8-gauss-points-required-for-second-order-hexahedral-finite-elements
"""
defaultquadrature(::Type{Edge2})   = GaussQuadrature{1}()
defaultquadrature(::Type{Edge3})   = GaussQuadrature{1}()
defaultquadrature(::Type{Tri3})    = GaussQuadrature{1}()
defaultquadrature(::Type{Tri6})    = GaussQuadrature{2}()
defaultquadrature(::Type{Quad4})   = GaussQuadrature{1}()
defaultquadrature(::Type{Quad8})   = GaussQuadrature{2}()
defaultquadrature(::Type{Tet4})    = GaussQuadrature{1}()
defaultquadrature(::Type{Tet10})   = GaussQuadrature{3}()
defaultquadrature(::Type{Hex8})    = GaussQuadrature{2}()
defaultquadrature(::Type{Hex20})   = GaussQuadrature{8}() # Why so high? 58 Points! Shouldn't 20 be enough?
defaultquadrature(::Type{Wedge6})  = GaussQuadrature{2}()
defaultquadrature(::Type{Wedge15}) = GaussQuadrature{8}() # Why so high? 46 Points! Shouldn't 15 be enough?

# Assume worst case for abstract types
defaultquadrature(::Type{<:AbstractShape1D}) = defaultquadrature(Edge3)
defaultquadrature(::Type{<:AbstractShape2D}) = defaultquadrature(Quad8)
defaultquadrature(::Type{<:AbstractShape3D}) = defaultquadrature(Hex20)

function defaultquadrature(sf::AbstractShapeFunctions)
    N = 2 * getorder(sf)
    return GaussQuadrature(N)
end

function defaultquadrature(sf::AbstractShapeFunctions, S::Type{<:AbstractShape})
    N = getorder(defaultquadrature(sf))
    M = getorder(defaultquadrature(S))
    return GaussQuadrature(max(N, M))
end
