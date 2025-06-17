export Lagrange
export dof_refcoordinates
export center_refcoordinates # This is actually has nothing to do with Lagrange...
export each_dof_refcoordinates
export edge2cellcoordinates
export face2cellcoordinates
export edge2elementcoordinates
export facet2elementcoordinates
export vertexdofs
export edgedofs
export facedofs
export facetdofs
export ndofs_geo

abstract type AbstractLagrange{N} <: AbstractScalarShapeFunctions{N} end

struct Lagrange{N} <: AbstractLagrange{N} end

Lagrange(order::Integer) = Lagrange{order}()

getorder(::Union{T, Type{T}}) where {T <: AbstractLagrange{N}} where {N} = N

sftype(::Type{<:Vertex})   = Lagrange{0}
sftype(::Type{<:Edge2})    = Lagrange{1}
sftype(::Type{<:Tri3})     = Lagrange{1}
sftype(::Type{<:Quad4})    = Lagrange{1}
sftype(::Type{<:Tet4})     = Lagrange{1}
sftype(::Type{<:Hex8})     = Lagrange{1}
sftype(::Type{<:Wedge6})   = Lagrange{1}
sftype(::Type{<:Edge3})    = Lagrange{2}
sftype(::Type{<:Tri6})     = Lagrange{2}
sftype(::Type{<:Quad8})    = Lagrange{2}
sftype(::Type{<:Tet10})    = Lagrange{2}
sftype(::Type{<:Hex20})    = Lagrange{2}
sftype(::Type{<:Wedge15})  = Lagrange{2}

shapetype(::Type{<:AbstractEdge},           ::Lagrange{1}) = Edge2
shapetype(::Type{<:AbstractEdge},           ::Lagrange{2}) = Edge3
shapetype(::Type{<:AbstractTriangle},       ::Lagrange{1}) = Tri3
shapetype(::Type{<:AbstractTriangle},       ::Lagrange{2}) = Tri6
shapetype(::Type{<:AbstractQuadrilateral},  ::Lagrange{1}) = Quad4
shapetype(::Type{<:AbstractQuadrilateral},  ::Lagrange{2}) = Quad8
shapetype(::Type{<:AbstractTetrahedron},    ::Lagrange{1}) = Tet4
shapetype(::Type{<:AbstractTetrahedron},    ::Lagrange{2}) = Tet10
shapetype(::Type{<:AbstractHexahedron},     ::Lagrange{1}) = Hex8
shapetype(::Type{<:AbstractHexahedron},     ::Lagrange{2}) = Hex20
shapetype(::Type{<:AbstractWedge},          ::Lagrange{1}) = Wedge6
shapetype(::Type{<:AbstractWedge},          ::Lagrange{2}) = Wedge15

ndofs(::Type{S}) where {S<:AbstractShape} = ndofs(S, sftype(S)())

edge2elementcoordinates( ::Type{T}, i::Integer, ξ::AbstractVector) where {T <: AbstractShape1D} = edge2edgecoordinates(T, i, ξ)
edge2elementcoordinates( ::Type{T}, i::Integer, ξ::AbstractVector) where {T <: AbstractShape2D} = edge2facecoordinates(T, i, ξ)
edge2elementcoordinates( ::Type{T}, i::Integer, ξ::AbstractVector) where {T <: AbstractShape3D} = edge2cellcoordinates(T, i, ξ)

facet2elementcoordinates(::Type{T}, i::Integer, ξ::AbstractVector) where {T <: AbstractShape1D} = vertex2edgecoordinates(T, i, ξ)
facet2elementcoordinates(::Type{T}, i::Integer, ξ::AbstractVector) where {T <: AbstractShape2D} = edge2facecoordinates(T, i, ξ)
facet2elementcoordinates(::Type{T}, i::Integer, ξ::AbstractVector) where {T <: AbstractShape3D} = face2cellcoordinates(T, i, ξ)

facetdofs(::Type{S}, sf::AbstractLagrange, i::Integer) where {S<:AbstractShape1D} = vertexdofs(S, sf, i)
facetdofs(::Type{S}, sf::AbstractLagrange, i::Integer) where {S<:AbstractShape2D} = edgedofs(S, sf, i)
facetdofs(::Type{S}, sf::AbstractLagrange, i::Integer) where {S<:AbstractShape3D} = facedofs(S, sf, i)

##############################
# Vertex 0D
##############################

ndofs(       ::Type{<:AbstractVertex}, ::Lagrange) = 1
ndofs_vertex(::Type{<:AbstractVertex}, ::Lagrange) = (1,)

function dof_refcoordinates(::Type{<:AbstractVertex}, ::Lagrange, i::Integer)
    i == 1 && return SVector{0, Float64}()
    throw(ArgumentError("No dof $i for Lagrange shape functions in AbstractVertex"))
end

function shapefunc_N(::Type{<:AbstractVertex}, ::Lagrange, i::Integer, ξ::AbstractVector)
    i == 1 && return 1.0
    throw(ArgumentError("No dof $i for Lagrange shape functions in AbstractVertex"))
end

function shapefunc_∇N(::Type{<:AbstractVertex}, ::Lagrange, i::Integer, ξ::AbstractVector)
    i == 1 && return SVector{0, Float64}()
    throw(ArgumentError("No dof $i for Lagrange shape functions in AbstractVertex"))
end

include("Lagrange1D.jl")
include("Lagrange2D.jl")
include("Lagrange3D.jl")
