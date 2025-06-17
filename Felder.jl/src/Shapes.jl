export AbstractShape, AbstractShape1D, AbstractShape2D, AbstractShape3D
export AbstractVertex, Vertex
export AbstractEdge, Edge2, Edge3
export AbstractFace
export AbstractTriangle, Tri3, Tri6
export AbstractQuadrilateral, Quad4, Quad8, Quad9
export AbstractCell
export AbstractTetrahedron, Tet4, Tet10
export AbstractHexahedron, Hex8, Hex20, Hex26
export AbstractWedge, Wedge6, Wedge15

export getid
export setid
export nnodes
export nvertices
export nedges
export nfaces
export ncells
export nfacets
export getnodes
export getvertexnodes
export getnormal
export gettangent
export getlocalvertices
export getlocaledges
export getlocalfaces
export getlocalvertex
export getlocaledge
export getlocalface
export localvertexindex
export localedgeindices
export localfaceindices
export localfacetindices
export localedgetype
export localedgetypes
export localfacetype
export localfacetypes
export localfacettype
export localfacettypes
export getlocalfacet
export flip
export getorientation
export getsubshape
export getsubshapes
export to_globalcoordinates
export center_coordinates
export parameterization_jacobian
export firstorder, firstordertype

##############################
# Shapes
##############################
# TODO: "local" indices should only be defined for the type.
# When dispatching on an instance of the type, we could return the "global" indices.

abstract type AbstractShape{Dim} end

# Type aliases
const AbstractShape0D = AbstractShape{0}
const AbstractShape1D = AbstractShape{1}
const AbstractShape2D = AbstractShape{2}
const AbstractShape3D = AbstractShape{3}

Base.ndims(::AbstractShape{Dim}) where {Dim} = Dim
Base.ndims(::Type{<:AbstractShape{Dim}}) where {Dim} = Dim

getid(s::AbstractShape) = s.id
getnodes(s::AbstractShape) = s.n

nnodes(::Union{T, Type{T}})     where {T <: AbstractShape} = nnodes(T)
nvertices(::Union{T, Type{T}})  where {T <: AbstractShape} = nvertices(T)
nedges(::Union{T, Type{T}})     where {T <: AbstractShape} = nedges(T)
nfaces(::Union{T, Type{T}})     where {T <: AbstractShape} = nfaces(T)
ncells(::Union{T, Type{T}})     where {T <: AbstractShape} = ncells(T)

nfacets(::Union{T, Type{T}})    where {T <: AbstractShape1D} = nvertices(T)
nfacets(::Union{T, Type{T}})    where {T <: AbstractShape2D} = nedges(T)
nfacets(::Union{T, Type{T}})    where {T <: AbstractShape3D} = nfaces(T)

# Default shape constructors
(::Type{T})(nodes::Int; id=0)             where {T <: AbstractShape} = T(tuple(nodes), id)
(::Type{T})(nodes::Tuple; id=0)           where {T <: AbstractShape} = T(nodes, id)
(::Type{T})(nodes::AbstractVector; id=0)  where {T <: AbstractShape} = T(nodes, id)
(::Type{T})(nodes...; id=0)               where {T <: AbstractShape} = T(nodes, id)

setid(s::T, id::Int) where {T<:AbstractShape} = T(s.n, id)

getvertexnodes(s::AbstractShape) = SVector{nvertices(s), Int}(s.n[localvertexindex(s, i)] for i in 1:nvertices(s))
getlocalvertex(s::AbstractShape, i::Integer) = Vertex(s.n[localvertexindex(s, i)])
getlocaledge(s::AbstractShape, i::Integer) = localedgetype(s, i)(s.n[localedgeindices(s, i)])
getlocalface(s::AbstractShape, i::Integer) = localfacetype(s, i)(s.n[localfaceindices(s, i)])

getlocalfacet(s::AbstractShape1D, i::Integer) = getlocalvertex(s, i)
getlocalfacet(s::AbstractShape2D, i::Integer) = getlocaledge(s, i)
getlocalfacet(s::AbstractShape3D, i::Integer) = getlocalface(s, i)

getlocalvertices(s::T) where {T <: AbstractShape} = ntuple(i -> getlocalvertex(s, i), Val(nvertices(T)))
getlocaledges(s::T) where {T <: AbstractShape}    = ntuple(i -> getlocaledge(s, i), Val(nedges(T)))
getlocalfaces(s::T) where {T <: AbstractShape}    = ntuple(i -> getlocalface(s, i), Val(nfaces(T)))

getlocalfacets(s::T) where {T <: AbstractShape1D} = getlocalvertices(s)
getlocalfacets(s::T) where {T <: AbstractShape2D} = getlocaledges(s)
getlocalfacets(s::T) where {T <: AbstractShape3D} = getlocalfaces(s)

localvertexindex(::T, i::Integer) where {T <: AbstractShape} = localvertexindex(T, i)
localedgeindices(::T, i::Integer) where {T <: AbstractShape} = localedgeindices(T, i)
localfaceindices(::T, i::Integer) where {T <: AbstractShape} = localfaceindices(T, i)

localedgetype(::Union{T, Type{T}}, i::Integer) where {T <: AbstractShape} = localedgetype(T, i)
localfacetype(::Union{T, Type{T}}, i::Integer) where {T <: AbstractShape} = localfacetype(T, i)

localedgetypes(::Union{T, Type{T}}) where {T <: AbstractShape} = ntuple(i -> localedgetype(T, i), Val(nedges(T)))
localfacetypes(::Union{T, Type{T}}) where {T <: AbstractShape} = ntuple(i -> localfacetype(T, i), Val(nfaces(T)))

localfacetindices(::Union{T, Type{T}}, i::Integer) where {T <: AbstractShape1D} = SA[localvertexindex(T, i)] # Put into SA for consistency with higher dimensions
localfacetindices(::Union{T, Type{T}}, i::Integer) where {T <: AbstractShape2D} = localedgeindices(T, i)
localfacetindices(::Union{T, Type{T}}, i::Integer) where {T <: AbstractShape3D} = localfaceindices(T, i)

localfacettype(::Union{T, Type{T}}, ::Integer) where {T <: AbstractShape1D} = Vertex
localfacettype(::Union{T, Type{T}}, i::Integer) where {T <: AbstractShape2D} = localedgetype(T, i)
localfacettype(::Union{T, Type{T}}, i::Integer) where {T <: AbstractShape3D} = localfacetype(T, i)

localfacettypes(::Union{T, Type{T}}) where {T <: AbstractShape2D} = localedgetypes(T)
localfacettypes(::Union{T, Type{T}}) where {T <: AbstractShape3D} = localfacetypes(T)

flip(s::T) where {T <: AbstractShape1D} = T(s.n[flipperm(T)], s.id)
flip(s::T) where {T <: AbstractShape2D} = T(s.n[flipperm(T)], s.id)
Base.circshift(s::T) where {T <: AbstractShape2D} = T(s.n[circshiftperm(T)], s.id)

nsubshapes(s::AbstractShape, ::Type{T}) where {T <: AbstractShape} = nsubshapes(typeof(s), T)
getsubshape(s::AbstractShape, ::Type{T}, i::Integer) where {T <: AbstractShape} = T(s.n[subshapeindices(typeof(s), T, i)])
getsubshapes(s::AbstractShape, ::Type{T}) where {T <: AbstractShape} = ntuple(i -> getsubshape(s, T, i), Val(nsubshapes(typeof(s), T)))

firstorder(s::T) where {T <: AbstractShape} = firstordertype(typeof(s))(s)

##############################
# Node Shapes
##############################

abstract type AbstractVertex <: AbstractShape{0} end

nnodes(::Type{<:AbstractVertex}) = 1
nvertices(::Type{<:AbstractVertex}) = 1

"""
Node defined by one node

    ○
    n1
"""
struct Vertex <: AbstractVertex
    n::SVector{1, Int}
    id::Int
end

function localvertexindex(::Type{Vertex}, i::Integer)
    i == 1 && return 1
    throw(ArgumentError("No local vertex $i for shape type Vertex."))
end

##############################
# Line Shapes
##############################

abstract type AbstractEdge <: AbstractShape{1} end

nvertices(::Type{<:AbstractEdge}) = 2
nedges(::Type{<:AbstractEdge}) = 1

"""
Linear line defined by 2 nodes

    ○────────────○
    n1           n2
"""
struct Edge2 <: AbstractEdge
    n::SVector{2, Int}
    id::Int
end

Edge2(n1::Int, n2::Int; id=0) = Edge2((n1, n2), id=id)

nnodes(::Type{Edge2}) = 2

function localvertexindex(::Type{Edge2}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    throw(ArgumentError("No local vertex $i for shape type Edge2."))
end

flipperm(::Type{Edge2}) = SA[2, 1]

Edge2(s::AbstractEdge) = Edge2(s.n[1], s.n[2], id=s.id)
firstordertype(::Type{<:AbstractEdge}) = Edge2

"""
Quadratic line defined by 3 nodes

    ○──────○──────○
    n1     n3     n2
"""
struct Edge3 <: AbstractEdge
    n::SVector{3, Int}
    id::Int
end

Edge3(n1::Int, n2::Int; id=0) = Edge3((n1, n2), id=id)

nnodes(::Type{Edge3}) = 3

function localvertexindex(::Type{Edge3}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    throw(ArgumentError("No local vertex $i for shape type Edge3."))
end

flipperm(::Type{Edge3}) = SA[2, 1 ,3]

##############################
# Face Shapes
##############################

abstract type AbstractFace <: AbstractShape{2} end

nfaces(::Type{<:AbstractFace}) = 1

# ------------- Triangle --------------

abstract type AbstractTriangle <: AbstractFace end

nvertices(::Type{<:AbstractTriangle}) = 3
nedges(::Type{<:AbstractTriangle}) = 3

"""
Linear triangle defined by 3 nodes

    n3
    ○
    ┊ ⋱
    ┊   ⋱
    ┊     ⋱
    ┊       ⋱
    ○┄┄┄┄┄┄┄┄┄○
    n1        n2
"""
struct Tri3 <: AbstractTriangle
    n::SVector{3, Int}
    id::Int
end

nnodes(::Type{Tri3}) = 3

function localvertexindex(::Type{Tri3}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    throw(ArgumentError("No local vertex $i for shape type Tri3."))
end

localedgetype(::Type{Tri3}, ::Integer) = Edge2

function localedgeindices(::Type{Tri3}, i::Integer)
    i == 1 && return SA[1, 2]
    i == 2 && return SA[2, 3]
    i == 3 && return SA[3, 1]
    throw(ArgumentError("No local edge $i for shape type Tri3."))
end

flipperm(::Type{Tri3})      = SA[3, 2, 1]
circshiftperm(::Type{Tri3}) = SA[2, 3, 1]

Tri3(s::AbstractTriangle) = Tri3(s.n[1], s.n[2], s.n[3], id=s.id)
firstordertype(::Type{<:AbstractTriangle}) = Tri3

"""
Quadratic triangle defined by 6 nodes

       n3
       ○
       ┊  ⋱
       ┊    ⋱
    n6 ○      ○ n5
       ┊        ⋱
       ┊          ⋱
       ○┄┄┄┄┄┄○┄┄┄┄┄┄○
       n1     n4     n2
"""
struct Tri6 <: AbstractTriangle
    n::SVector{6, Int}
    id::Int
end

nnodes(::Type{Tri6}) = 6

function localvertexindex(::Type{Tri6}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    throw(ArgumentError("No local vertex $i for shape type Tri6."))
end

localedgetype(::Type{Tri6}, ::Integer) = Edge3

# Defined such that `getnormal(edge)` points outwards in 2D
function localedgeindices(::Type{Tri6}, i::Integer)
    i == 1 && return SA[1, 2, 4]
    i == 2 && return SA[2, 3, 5]
    i == 3 && return SA[3, 1, 6]
    throw(ArgumentError("No local edge $i for shape type Tri6."))
end

flipperm(::Type{Tri6})      = SA[3, 2, 1, 5, 4, 6]
circshiftperm(::Type{Tri6}) = SA[2, 3, 1, 5, 6, 4]

# ------------- Quadrilateral --------------

abstract type AbstractQuadrilateral <: AbstractFace end

nvertices(::Type{<:AbstractQuadrilateral}) = 4
nedges(::Type{<:AbstractQuadrilateral}) = 4

"""
Linear quadrilateral defined by 4 nodes

    n4            n3
    ○─────────────○
    │             │
    │             │
    │             │
    ○─────────────○
    n1            n2
"""
struct Quad4 <: AbstractQuadrilateral
    n::SVector{4, Int}
    id::Int
end

nnodes(::Type{Quad4}) = 4

function localvertexindex(::Type{Quad4}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    throw(ArgumentError("No local vertex $i for shape type Quad4."))
end

localedgetype(::Type{Quad4}, ::Integer) = Edge2

# Defined such that `getnormal(edge)` points outwards in 2D
function localedgeindices(::Type{Quad4}, i::Integer)
    i == 1 && return SA[1, 2]
    i == 2 && return SA[2, 3]
    i == 3 && return SA[3, 4]
    i == 4 && return SA[4, 1]
    throw(ArgumentError("No local edge $i for shape type Quad4."))
end

flipperm(::Type{Quad4})      = SA[4, 3, 2, 1]
circshiftperm(::Type{Quad4}) = SA[2, 3, 4, 1]

Quad4(s::AbstractQuadrilateral) = Quad4(s.n[1], s.n[2], s.n[3], s.n[4], id=s.id)
firstordertype(::Type{<:AbstractQuadrilateral}) = Quad4

"""
Quadratic quadrilateral defined by 8 nodes

    n4     n7     n3
    ○──────○──────○
    │             │
    ○ n8          ○ n6
    │             │
    ○──────○──────○
    n1     n5     n2
"""
struct Quad8 <: AbstractQuadrilateral
    n::SVector{8, Int}
    id::Int
end

nnodes(::Type{Quad8}) = 8

function localvertexindex(::Type{Quad8}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    throw(ArgumentError("No local vertex $i for shape type Quad8."))
end

localedgetype(::Type{Quad8}, ::Integer) = Edge3

# Defined such that `getnormal(edge)` points outwards in 2D
function localedgeindices(::Type{Quad8}, i::Integer)
    i == 1 && return SA[1, 2, 5]
    i == 2 && return SA[2, 3, 6]
    i == 3 && return SA[3, 4, 7]
    i == 4 && return SA[4, 1, 8]
    throw(ArgumentError("No local edge $i for shape type Quad8."))
end

flipperm(::Type{Quad8})      = SA[4, 3, 2, 1, 7, 6, 5, 8]
circshiftperm(::Type{Quad8}) = SA[2, 3, 4, 1, 6, 7, 8, 5]

"""
Quadratic quadrilateral defined by 9 nodes.

Used for mesh generation of structured grids, because
it can be easily divided into two Tri6 shapes.

    n4     n7     n3
    ○──────○──────○
    │             │
    ○ n8   o n9   ○ n6
    │             │
    ○──────○──────○
    n1     n5     n2
"""
struct Quad9 <: AbstractQuadrilateral
    n::SVector{9, Int}
    id::Int
end

nnodes(::Type{Quad9}) = 9

function localvertexindex(::Type{Quad9}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    throw(ArgumentError("No local vertex $i for shape type Quad9."))
end

localedgetype(::Type{Quad9}, ::Integer) = Edge3

# Defined such that `getnormal(edge)` points outwards in 2D
function localedgeindices(::Type{Quad9}, i::Integer)
    i == 1 && return SA[1, 2, 5]
    i == 2 && return SA[2, 3, 6]
    i == 3 && return SA[3, 4, 7]
    i == 4 && return SA[4, 1, 8]
    throw(ArgumentError("No local edge $i for shape type Quad9."))
end

flipperm(::Type{Quad9})      = SA[4, 3, 2, 1, 7, 6, 5, 8, 9]
circshiftperm(::Type{Quad9}) = SA[2, 3, 4, 1, 6, 7, 8, 5, 9]

##############################
# Volume Shapes
##############################

abstract type AbstractCell <: AbstractShape{3} end

# ------------- Tetrahedron --------------

abstract type AbstractTetrahedron <: AbstractCell end

nvertices(::Type{<:AbstractTetrahedron}) = 4
nedges(::Type{<:AbstractTetrahedron}) = 6
nfaces(::Type{<:AbstractTetrahedron}) = 4

"""
Linear Tetrahedron defined by 4 nodes

    n4
    ○
    │╲⋱
    │ ╲ ⋱
    │  ╲  ⋱
    │n3 ○   ⋱
    │  /  ⋱   ⋱
    │ /      ⋱  ⋱
    │/          ⋱ ⋱
    ○───────────────○
    n1              n2

"""
struct Tet4 <: AbstractTetrahedron
    n::SVector{4, Int}
    id::Int
end

nnodes(::Type{Tet4}) = 4

function localvertexindex(::Type{Tet4}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    throw(ArgumentError("No local vertex $i for shape type Tet4."))
end

localedgetype(::Type{Tet4}, ::Integer) = Edge2

function localedgeindices(::Type{Tet4}, i::Integer)
    i == 1 && return SA[1, 2]
    i == 2 && return SA[2, 3]
    i == 3 && return SA[3, 1]
    i == 4 && return SA[1, 4]
    i == 5 && return SA[2, 4]
    i == 6 && return SA[3, 4]
    throw(ArgumentError("No local edge $i for shape type Tet4."))
end

localfacetype(::Type{Tet4}, ::Integer) = Tri3

# Defined such that `getnormal(face)` points outwards
function localfaceindices(::Type{Tet4}, i::Integer)
    i == 1 && return SA[1, 2, 4]
    i == 2 && return SA[2, 3, 4]
    i == 3 && return SA[3, 1, 4]
    i == 4 && return SA[1, 3, 2]
    throw(ArgumentError("No local face $i for shape type Tet4."))
end

Tet4(s::AbstractTetrahedron) = Tet4(s.n[1], s.n[2], s.n[3], s.n[4], id=s.id)
firstordertype(::Type{<:AbstractTetrahedron}) = Tet4

"""
Quadratic Tetrahedron defined by 10 nodes

       n4
       ○
       │╲⋱
    n10│ ○ ⋱
       │  ╲  ⋱
    n8 ○   ○   ○ n9
       │  /n3 ⋱  ⋱
     n7│ ○   n6 ○  ⋱
       │/          ⋱ ⋱
       ○───────○───────○
       n1      n5      n2

"""
struct Tet10 <: AbstractTetrahedron
    n::SVector{10, Int}
    id::Int
end

nnodes(::Type{Tet10}) = 10

function localvertexindex(::Type{Tet10}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    throw(ArgumentError("No local vertex $i for shape type Tet10."))
end

localedgetype(::Type{Tet10}, ::Integer) = Edge3

function localedgeindices(::Type{Tet10}, i::Integer)
    i == 1 && return SA[1, 2, 5]
    i == 2 && return SA[2, 3, 6]
    i == 3 && return SA[3, 1, 7]
    i == 4 && return SA[1, 4, 8]
    i == 5 && return SA[2, 4, 9]
    i == 6 && return SA[3, 4, 10]
    throw(ArgumentError("No local edge $i for shape type Tet10."))
end

localfacetype(::Type{Tet10}, ::Integer) = Tri6

# Defined such that `getnormal(face)` points outwards
function localfaceindices(::Type{Tet10}, i::Integer)
    i == 1 && return SA[1, 2, 4, 5,  9,  8]
    i == 2 && return SA[2, 3, 4, 6, 10,  9]
    i == 3 && return SA[3, 1, 4, 7,  8, 10]
    i == 4 && return SA[1, 3, 2, 7,  6,  5]
    throw(ArgumentError("No local face $i for shape type Tet10."))
end

# ------------- Hexahedron --------------

abstract type AbstractHexahedron <: AbstractCell end

nvertices(::Type{<:AbstractHexahedron}) = 8
nedges(::Type{<:AbstractHexahedron}) = 12
nfaces(::Type{<:AbstractHexahedron}) = 6

"""
Linear Hexahedron defined by 8 nodes

           n8            n7
           ○─────────────○
          /│            /│
         / │           / │
        /  │          /  │
    n5 ○─────────────○ n6│
       │   │         │   │
       │n4 ○─────────│───○ n3
       │  /          │  /
       │ /           │ /
       │/            │/
       ○─────────────○
       n1            n2
"""
struct Hex8 <: AbstractHexahedron
    n::SVector{8, Int}
    id::Int
end

nnodes(::Type{Hex8}) = 8

function localvertexindex(::Type{Hex8}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    i == 5 && return 5
    i == 6 && return 6
    i == 7 && return 7
    i == 8 && return 8
    throw(ArgumentError("No local vertex $i for shape type Hex8."))
end

localedgetype(::Type{Hex8}, ::Integer) = Edge2

function localedgeindices(::Type{Hex8}, i::Integer)
    i ==  1 && return SA[1, 2]
    i ==  2 && return SA[2, 3]
    i ==  3 && return SA[3, 4]
    i ==  4 && return SA[4, 1]
    i ==  5 && return SA[5, 6]
    i ==  6 && return SA[6, 7]
    i ==  7 && return SA[7, 8]
    i ==  8 && return SA[8, 5]
    i ==  9 && return SA[1, 5]
    i == 10 && return SA[2, 6]
    i == 11 && return SA[3, 7]
    i == 12 && return SA[4, 8]
    throw(ArgumentError("No local edge $i for shape type Hex8."))
end

localfacetype(::Type{Hex8}, ::Integer) = Quad4

# Defined such that `getnormal(face)` points outwards
function localfaceindices(::Type{Hex8}, i::Integer)
    i == 1 && return SA[1, 2, 6, 5]
    i == 2 && return SA[2, 3, 7, 6]
    i == 3 && return SA[3, 4, 8, 7]
    i == 4 && return SA[4, 1, 5, 8]
    i == 5 && return SA[1, 4, 3, 2]
    i == 6 && return SA[5, 6, 7, 8]
    throw(ArgumentError("No local face $i for shape type Hex8."))
end

Hex8(s::AbstractHexahedron) = Hex8(s.n[1], s.n[2], s.n[3], s.n[4],
    s.n[5], s.n[6], s.n[7], s.n[8], id=s.id)
firstordertype(::Type{<:AbstractHexahedron}) = Hex8

"""
Quadratic Hexahedron defined by 20 nodes

            n8      n15    n7
            ○──────○──────○
           /│            /│
      n16 ○ │       n14 ○ │
         /  ○ n20      /  ○ n19
     n5 ○──────○──────○ n6│
        │   │ n13 n11 │   │
        │n4 ○──────○──│───○ n3
    n17 ○  /       n18○  /
        │ ○ n12       │ ○ n10
        │/            │/
        ○──────○──────○
        n1      n9     n2
"""
struct Hex20 <: AbstractHexahedron
    n::SVector{20, Int}
    id::Int
end

nnodes(::Type{Hex20}) = 20

function localvertexindex(::Type{Hex20}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    i == 5 && return 5
    i == 6 && return 6
    i == 7 && return 7
    i == 8 && return 8
    throw(ArgumentError("No local vertex $i for shape type Hex20."))
end

localedgetype(::Type{Hex20}, ::Integer) = Edge3

function localedgeindices(::Type{Hex20}, i::Integer)
    i ==  1 && return SA[1, 2,  9]
    i ==  2 && return SA[2, 3, 10]
    i ==  3 && return SA[3, 4, 11]
    i ==  4 && return SA[4, 1, 12]
    i ==  5 && return SA[5, 6, 13]
    i ==  6 && return SA[6, 7, 14]
    i ==  7 && return SA[7, 8, 15]
    i ==  8 && return SA[8, 5, 16]
    i ==  9 && return SA[1, 5, 17]
    i == 10 && return SA[2, 6, 18]
    i == 11 && return SA[3, 7, 19]
    i == 12 && return SA[4, 8, 20]
    throw(ArgumentError("No local edge $i for shape type Hex20."))
end

localfacetype(::Type{Hex20}, ::Integer) = Quad8

# Defined such that `getnormal(face)` points outwards
function localfaceindices(::Type{Hex20}, i::Integer)
    i == 1 && return SA[1, 2, 6, 5,  9, 18, 13, 17]
    i == 2 && return SA[2, 3, 7, 6, 10, 19, 14, 18]
    i == 3 && return SA[3, 4, 8, 7, 11, 20, 15, 19]
    i == 4 && return SA[4, 1, 5, 8, 12, 17, 16, 20]
    i == 5 && return SA[1, 4, 3, 2, 12, 11, 10,  9]
    i == 6 && return SA[5, 6, 7, 8, 13, 14, 15, 16]
    throw(ArgumentError("No local face $i for shape type Hex20."))
end

"""
Quadratic Hexahedron defined by 26 nodes
(8 vertex nodes, 12 edge nodes and 6 face nodes).

Used for mesh generation of structured 3D grids, because
it can be easily divided into five Tet10 shapes or
or two Wedge15 shapes.

            n8      n15    n7
            ○────────○────────○
           /│     n26        /│
      n16 ○ │      o    n14 ○ │
         /  ○ n20    o n23 /  ○ n19
     n5 ○────────○────────○ n6│
      n24 o │ n13     n11 │ o n22
        │n4 ○────────○────│───○ n3
    n17 ○  / n21 o    n18 ○  /
      n12 ○        o      │ ○ n10
        │/        n25     │/
        ○────────○────────○
        n1      n9     n2
"""
struct Hex26 <: AbstractHexahedron
    n::SVector{26, Int}
    id::Int
end

nnodes(::Type{Hex26}) = 26

function localvertexindex(::Type{Hex26}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    i == 5 && return 5
    i == 6 && return 6
    i == 7 && return 7
    i == 8 && return 8
    throw(ArgumentError("No local vertex $i for shape type Hex26."))
end

localedgetype(::Type{Hex26}, ::Integer) = Edge3

function localedgeindices(::Type{Hex26}, i::Integer)
    i ==  1 && return SA[1, 2,  9]
    i ==  2 && return SA[2, 3, 10]
    i ==  3 && return SA[3, 4, 11]
    i ==  4 && return SA[4, 1, 12]
    i ==  5 && return SA[5, 6, 13]
    i ==  6 && return SA[6, 7, 14]
    i ==  7 && return SA[7, 8, 15]
    i ==  8 && return SA[8, 5, 16]
    i ==  9 && return SA[1, 5, 17]
    i == 10 && return SA[2, 6, 18]
    i == 11 && return SA[3, 7, 19]
    i == 12 && return SA[4, 8, 20]
    throw(ArgumentError("No local edge $i for shape type Hex26."))
end

localfacetype(::Type{Hex26}, ::Integer) = Quad9

# Defined such that `getnormal(face)` points outwards
function localfaceindices(::Type{Hex26}, i::Integer)
    i == 1 && return SA[1, 2, 6, 5,  9, 18, 13, 17, 21]
    i == 2 && return SA[2, 3, 7, 6, 10, 19, 14, 18, 22]
    i == 3 && return SA[3, 4, 8, 7, 11, 20, 15, 19, 23]
    i == 4 && return SA[4, 1, 5, 8, 12, 17, 16, 20, 24]
    i == 5 && return SA[1, 4, 3, 2, 12, 11, 10,  9, 25]
    i == 6 && return SA[5, 6, 7, 8, 13, 14, 15, 16, 26]
    throw(ArgumentError("No local face $i for shape type Hex26."))
end

# ------------- Wedge --------------

abstract type AbstractWedge <: AbstractCell end

nvertices(::Type{<:AbstractWedge}) = 6
nedges(::Type{<:AbstractWedge}) = 9
nfaces(::Type{<:AbstractWedge}) = 5

"""
Linear Wedge defined by 6 points

        p6 ○
          /│  ⋱
         / │     ⋱
        /  │        ⋱
    p4 ○───────────────○ p5
       │   │           │
       │   │           │
       │   │           │
       │   ○ p3        │
       │  /   ⋱        │
       │ /       ⋱     │
       │/           ⋱  │
       ○───────────────○
       p1              p2

"""
struct Wedge6 <: AbstractWedge
    n::SVector{6, Int}
    id::Int
end

nnodes(::Type{<:Wedge6}) = 6

function localvertexindex(::Type{Wedge6}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    i == 5 && return 5
    i == 6 && return 6
    throw(ArgumentError("No local vertex $i for shape type Wedge6."))
end

localedgetype(::Type{Wedge6}, ::Integer) = Edge2

function localedgeindices(::Type{Wedge6}, i::Integer)
    i == 1 && return SA[1, 2]
    i == 2 && return SA[2, 3]
    i == 3 && return SA[3, 1]
    i == 4 && return SA[4, 5]
    i == 5 && return SA[5, 6]
    i == 6 && return SA[6, 4]
    i == 7 && return SA[1, 4]
    i == 8 && return SA[2, 5]
    i == 9 && return SA[3, 6]
    throw(ArgumentError("No local edge $i for shape type Wedge6."))
end

function localfacetype(::Type{Wedge6}, i::Integer)
    i == 1 && return Quad4
    i == 2 && return Quad4
    i == 3 && return Quad4
    i == 4 && return Tri3
    i == 5 && return Tri3
    throw(ArgumentError("No local face $i for shape type Wedge6."))
end

# Defined such that `getnormal(face)` points outwards
# Type instability!
function localfaceindices(::Type{Wedge6}, i::Integer)
    i == 1 && return SA[1, 2, 5, 4]
    i == 2 && return SA[2, 3, 6, 5]
    i == 3 && return SA[3, 1, 4, 6]
    i == 4 && return SA[1, 3, 2]
    i == 5 && return SA[4, 5, 6]
    throw(ArgumentError("No local face $i for shape type Wedge6."))
end

Wedge6(s::AbstractWedge) = Wedge6(s.n[1], s.n[2], s.n[3], s.n[4], s.n[5], s.n[6], id=s.id)
firstordertype(::Type{<:AbstractWedge}) = Wedge6

"""
Quadratic Wedge defined by 15 points

        p6 ○
          /│  ⋱
     p12 ○ │     ○ p11
        /  │        ⋱
    p4 ○───o───○───────○ p5
       │p15│   p10     │
       │   │           │
       │   │           │
   p13 ○   ○ p3        ○ p14
       │  /   ⋱        │
       │ ○ p9    ○ p8  │
       │/           ⋱  │
       ○───────○───────○
       p1      p7      p2
"""
struct Wedge15 <: AbstractWedge
    n::SVector{15, Int}
    id::Int
end

nnodes(::Type{<:Wedge15}) = 15

function localvertexindex(::Type{Wedge15}, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    i == 5 && return 5
    i == 6 && return 6
    throw(ArgumentError("No local vertex $i for shape type Wedge15."))
end

localedgetype(::Type{Wedge15}, ::Integer) = Edge3

function localedgeindices(::Type{Wedge15}, i::Integer)
    i == 1 && return SA[1, 2, 7]
    i == 2 && return SA[2, 3, 8]
    i == 3 && return SA[3, 1, 9]
    i == 4 && return SA[4, 5, 10]
    i == 5 && return SA[5, 6, 11]
    i == 6 && return SA[6, 4, 12]
    i == 7 && return SA[1, 4, 13]
    i == 8 && return SA[2, 5, 14]
    i == 9 && return SA[3, 6, 15]
    throw(ArgumentError("No local edge $i for shape type Wedge15."))
end

function localfacetype(::Type{Wedge15}, i::Integer)
    i == 1 && return Quad8
    i == 2 && return Quad8
    i == 3 && return Quad8
    i == 4 && return Tri6
    i == 5 && return Tri6
    throw(ArgumentError("No local face $i for shape type Wedge15."))
end

# Defined such that `getnormal(face)` points outwards
function localfaceindices(::Type{Wedge15}, i::Integer)
    i == 1 && return SA[1, 2, 5, 4, 7, 14, 10, 13]
    i == 2 && return SA[2, 3, 6, 5, 8, 15, 11, 14]
    i == 3 && return SA[3, 1, 4, 6, 9, 13, 12, 15]
    i == 4 && return SA[1, 3, 2, 9, 8, 7]
    i == 5 && return SA[4, 5, 6, 10, 11, 12]
    throw(ArgumentError("No local face $i for shape type Wedge15."))
end

##############################
# Parameterization
##############################

"""
    to_globalcoordinates(ξ, shapetype, shapecoordinates)

Convert the reference coordinates `ξ` on the shape of type `shapetype` defined by
`shapecoordinates` to global coordinates.

WARNING: It is not checked if reference coordinates ξ actually lie inside the
reference shape!

# Examples
```julia-repl
julia> shapecoordinates = [
           SA[1.0, 1.0],
           SA[2.0, 1.0],
           SA[1.0, 2.0],
           SA[1.5, 1.0],
           SA[1.7, 1.7],
           SA[1.28, 1.5]
       ];

julia> to_globalcoordinates([0, 0], Tri6, shapecoordinates)
2-element SVector{2, Float64} with indices SOneTo(2):
 1.0
 1.0

julia> to_globalcoordinates([1, 0], Tri6, shapecoordinates)
2-element SVector{2, Float64} with indices SOneTo(2):
 2.0
 1.0

julia> to_globalcoordinates([0.5, 0.5], Tri6, shapecoordinates)
2-element SVector{2, Float64} with indices SOneTo(2):
 1.7
 1.7
```
"""
function to_globalcoordinates(ξ, ::Type{T}, shapecoordinates) where {T <: AbstractShape}
    shapefunc = sftype(T)()
    @assert ndofs(T, shapefunc) == length(shapecoordinates)
    x = zero(eltype(shapecoordinates))
    for (i, Ni) in enumerate(eachshapefunc_N(T, shapefunc, ξ))
        @inbounds x += Ni * shapecoordinates[i]
    end
    return x
end

"""
    center_coordinates(shapetype, shapecoordinates)
    center_coordinates(coordinates, shapes::Vector{<:AbstractShape})

Return the global coordinates at the shape center as defined by
`center_refcoordinates(shapetype.)`.

Note: This is probably not the center of mass in general.
"""
function center_coordinates(::Type{T}, shapecoordinates) where {T <: AbstractShape}
    to_globalcoordinates(center_refcoordinates(T), T, shapecoordinates)
end

function center_coordinates(coordinates, shapes::Vector{T}) where {T <: AbstractShape}
    [center_coordinates(typeof(shape), view(coordinates, shape.n)) for shape in shapes]
end

"""
    parameterization_jacobian(ξ, shapetype, shapecoordinates)

Get the Jacobian of the parametrization from reference to global coordinates
at reference coordinates `ξ` in a shape of type `shapetype`.

For example, if the parameterization of a surface with reference coordinates ξ = [s, t]
in 3D space is r(s, t) = [r1(s, t), r2(s, t)], the 3x2 Jacobian

    J(r) = [∂r/∂s ∂r/∂t]

is returned.

WARNING: It is not checked if reference coordinates ξ actually lie inside the
reference shape!

TODO: The Jacobian is constant for shapes of geometric order 1. Simplify calculation?

# Examples
```julia-repl
julia> shapecoordinates = [
           SA[1.0, 1.0],
           SA[2.0, 1.0],
           SA[1.0, 2.0],
           SA[1.5, 1.0],
           SA[1.7, 1.7],
           SA[1.28, 1.5]
       ];

julia> parameterization_jacobian([0, 0], Tri6, shapecoordinates)
2×2 SMatrix{2, 2, Float64, 4} with indices SOneTo(2)×SOneTo(2):
 1.0  1.12
 0.0  1.0

julia> parameterization_jacobian([0.5, 0.5], Tri6, shapecoordinates)
2×2 SMatrix{2, 2, Float64, 4} with indices SOneTo(2)×SOneTo(2):
 0.84  -0.16
 0.4    1.4

julia> shapecoordinates = [
           SA[1.0, 1.0, 1.0],
           SA[2.0, 1.0, 1.0],
           SA[1.0, 2.0, 1.0],
           SA[1.5, 1.0, 1.0],
           SA[1.7, 1.7, 1.0],
           SA[1.28, 1.5, 1.0]
       ];

julia> parameterization_jacobian([0, 0], Tri6, shapecoordinates)
3×2 SMatrix{3, 2, Float64, 6} with indices SOneTo(3)×SOneTo(2):
 1.0  1.12
 0.0  1.0
 0.0  0.0
```
"""
function parameterization_jacobian(ξ, ::Type{T}, shapecoordinates) where {T <: AbstractShape}
    @assert ndofs(T) == length(shapecoordinates)
    @assert ndims(T) == length(ξ)
    N = length(eltype(shapecoordinates))
    M = ndims(T)
    J = zero(SMatrix{N, M, Float64})
    for (xi, dNdξi) in zip(shapecoordinates, eachshapefunc_∇N(T, sftype(T)(), ξ))
        @inbounds J += xi .* dNdξi'
    end
    return J
end

"""
    parameterization_normal(ξ, shapetype0D, shapecoordinates1D)
    parameterization_normal(shapetype0D, jacobian1x0)

Return the normal vector of the parameterization of a vertex (0D shape) in 1D space
at reference coordinates `ξ`, which here defined as a 1-element NaN vector.

# Note: This is not a well-defined concept in 1D space and is only implemented
here for consistency with other shapes in higher dimensions.
"""
function parameterization_normal(ξ, ::Type{T}, shapecoordinates) where {T <: AbstractShape0D}
    J = parameterization_jacobian(ξ, T, shapecoordinates)
    return parameterization_normal(T, J)
end

function parameterization_normal(::Type{<:AbstractShape0D}, J::AbstractMatrix)
    @assert size(J) == (1, 0)
    return SA[1.0]
end

"""
    parameterization_normal(ξ, shapetype1D, shapecoordinates2D)
    parameterization_normal(shapetype1D, jacobian2x1)

Return the normal vector of the parameterization of a 1D shape in 2D space
at reference coordinates `ξ`.

NOTE: The returned vector does NOT have unit length so it can be used as the
differential line element in line integrals??? TODO: Check this! What is
the orientation of the normal on the line in 2D space? Left, right?

# Note: A normal vector on a 1D line in 3D space is not well-defined.

# References
https://en.wikipedia.org/wiki/Line_integral
"""
function parameterization_normal(ξ, ::Type{T}, shapecoordinates,) where {T <: AbstractShape1D}
    J = parameterization_jacobian(ξ, T, shapecoordinates)
    return parameterization_normal(T, J)
end

function parameterization_normal(::Type{<:AbstractShape1D}, J::AbstractMatrix)
    @assert size(J) == (2, 1)
    return SA[J[2, 1], -J[1, 1]]
end

"""
    parameterization_normal(ξ, shapetype1D, shapecoordinates1D, ivertex)
    parameterization_normal(shapetype1D, jacobian1x1, ivertex)

A normal vector on a 1D line in 1D space is not defined, so this function
returns the sign ("-1" or "1") of the "direction" of the vertex `ivertex` on the
line element as a SVector{1}. "1" means it is the "right" vertex,
"-1" means it is the "left" vertex.

The method has the name `parameterization_normal` to be consistent with other
shapes in higher dimensions.

WARNING: It is not checked if reference coordinates `ξ` actualle lie on edge `ivertex`.

# References
https://en.wikipedia.org/wiki/Line_integral
"""
function parameterization_normal(ξ, ::Type{T}, shapecoordinates, ivertex) where {T <: AbstractShape1D}
    J = parameterization_jacobian(ξ, T, shapecoordinates)
    return parameterization_normal(T, J, ivertex)
end

function parameterization_normal(::Type{<:AbstractEdge}, J::AbstractMatrix, ivertex::Integer)
    @assert size(J) == (1, 1)
    ivertex == 1 && return SA[-one(eltype(J))]
    ivertex == 2 && return SA[ one(eltype(J))]
    throw(ArgumentError("no normal defined for vertex $ivertex in AbstractEdge"))
end

"""
    parameterization_normal(ξ, shapetype2D, shapecoordinates3D)
    parameterization_normal(shapetype2D, jacobian3x2)

Return the normal vector of the parameterization of a 2D shape in 3D space
at reference coordinates `ξ`.

NOTE: The returned vector does NOT have unit length so it can be used as the
differential surface element in surface integrals.

# References
https://en.wikipedia.org/wiki/Surface_integral
"""
function parameterization_normal(ξ, ::Type{T}, shapecoordinates) where {T <: AbstractShape2D}
    J = parameterization_jacobian(ξ, T, shapecoordinates)
    return parameterization_normal(T, J)
end

function parameterization_normal(::Type{<:AbstractShape2D}, J::AbstractMatrix)
    @assert size(J) == (3, 2)
    return cross(J[:, 1], J[:, 2])
end

"""
    parameterization_normal(ξ, shapetype2D, shapecoordinates2D, iedge)
    parameterization_normal(shapetype2D, jacobian2x2, iedge)

Return the (outwards) normal vector of edge `iedge` of the parameterization of
a 2D shape in 2D space at reference coordinates `ξ`. `ξ` must describe a point on
the local edge `iedge`.

WARNING: It is not checked if reference coordinates `ξ` actualle lie on edge `iedge`.

NOTE: The returned vector does NOT have unit length so it can be used as the
differential surface element in surface integrals.

# References
https://en.wikipedia.org/wiki/Surface_integral
"""
function parameterization_normal(ξ, ::Type{T}, shapecoordinates, iedge) where {T <: AbstractShape2D}
    J = parameterization_jacobian(ξ, T, shapecoordinates)
    return parameterization_normal(T, J, iedge)
end

function parameterization_normal(::Type{<:AbstractTriangle}, J::AbstractMatrix, iedge::Integer)
    @assert size(J) == (2, 2)
    iedge == 1 && return SA[J[2, 1], -J[1, 1]]
    iedge == 2 && return SA[-J[2, 1] + J[2, 2], J[1, 1] - J[1, 2]]
    iedge == 3 && return SA[-J[2, 2], J[1, 2]]
    throw(ArgumentError("no normal defined for edge $iedge in AbstactTriangle"))
end

function parameterization_normal(::Type{<:AbstractQuadrilateral}, J::AbstractMatrix, iedge::Integer)
    @assert size(J) == (2, 2)
    iedge == 1 && return SA[J[2, 1], -J[1, 1]]
    iedge == 2 && return SA[J[2, 2], -J[1, 2]]
    iedge == 3 && return SA[-J[2, 1], J[1, 1]]
    iedge == 4 && return SA[-J[2, 2], J[1, 2]]
    throw(ArgumentError("no normal defined for edge $iedge in AbstractQuadrilateral"))
end

"""
    parameterization_normal(ξ, shapetype3D, shapecoordinates3D, iface)
    parameterization_normal(shapetype3D, jacobian3x3, iface)

Return the (outwards) normal vector of face `iface` of the parameterization of
a 3D shape in 3D space at reference coordinates `ξ`. `ξ` must describe a point on
the local face `iface`.

WARNING: It is not checked if reference coordinates `ξ` actualle lie on face `iface`.

NOTE: The returned vector does NOT have unit length so it can be used as the
differential surface element in surface integrals.

# References
https://en.wikipedia.org/wiki/Surface_integral
"""
function parameterization_normal(ξ, ::Type{T}, shapecoordinates, iface) where {T <: AbstractShape3D}
    J = parameterization_jacobian(ξ, T, shapecoordinates)
    return parameterization_normal(T, J, iface)
end

function parameterization_normal(::Type{<:AbstractTetrahedron}, J::AbstractMatrix, iface::Integer)
    @assert size(J) == (3, 3)
    iface == 1 && return cross(J[:, 1], J[:, 3])
    iface == 2 && return cross(-J[:, 1] + J[:, 2], -J[:, 1] + J[:, 3])
    iface == 3 && return cross(J[:, 3], J[:, 2])
    iface == 4 && return cross(J[:, 2], J[:, 1])
    throw(ArgumentError("no normal defined for face $iface in AbstractTetrahedron"))
end

function parameterization_normal(::Type{<:AbstractHexahedron}, J::AbstractMatrix, iface::Integer)
    @assert size(J) == (3, 3)
    iface == 1 && return cross(J[:, 1], J[:, 3])
    iface == 2 && return cross(J[:, 2], J[:, 3])
    iface == 3 && return cross(J[:, 3], J[:, 1])
    iface == 4 && return cross(J[:, 3], J[:, 2])
    iface == 5 && return cross(J[:, 2], J[:, 1])
    iface == 6 && return cross(J[:, 1], J[:, 2])
    throw(ArgumentError("no normal defined for face $iface in AbstractHexahedron"))
end

function parameterization_normal(::Type{<:AbstractWedge}, J::AbstractMatrix, iface::Integer)
    @assert size(J) == (3, 3)
    iface == 1 && return cross(J[:, 1], J[:, 3])
    iface == 2 && return cross(-J[:, 1] + J[:, 2], J[:, 3])
    iface == 3 && return cross(J[:, 3], J[:, 2])
    iface == 4 && return cross(J[:, 2], J[:, 1])
    iface == 5 && return cross(J[:, 1], J[:, 2])
    throw(ArgumentError("no normal defined for face $iface in AbstractWedge"))
end

"""
    getnormal(ξ, shapetype2D, shapecoordinates3D)

Return the normal vector (unit length) of the parameterization of a 2D shape
in 3D space at reference coordinates `ξ`.

TODO: Shapes of first geometric order have a constant normal and do not need ξ.
"""
function getnormal(ξ, ::Type{T}, shapecoordinates) where {T <: AbstractShape2D}
    n = parameterization_normal(ξ, T, shapecoordinates)
    return n / norm(n)
end

"""
    getnormal(ξ, shapetype2D, shapecoordinates2D, iedge)

Return the (outward) normal vector (unit length) of edge `iedge` of the parameterization
of a 2D shape in 2D space at reference coordinates `ξ`. `ξ` must describe a point on
the local edge `iedge`.

WARNING: It is not checked if reference coordinates `ξ` actualle lie on face `iedge`.

TODO: Shapes of first geometric order have a constant normal and do not need ξ.
"""
function getnormal(ξ, ::Type{T}, shapecoordinates, iedge) where {T <: AbstractShape2D}
    n = parameterization_normal(ξ, T, shapecoordinates, iedge)
    return n / norm(n)
end

"""
    getnormal(ξ, shapetype3D, shapecoordinates3D, iface)

Return the (outward) normal vector (unit length) of face `iface` of the parameterization
of a 3D shape in 3D space at reference coordinates `ξ`. `ξ` must describe a point on
the local face `iface`.

WARNING: It is not checked if reference coordinates `ξ` actualle lie on face `iface`.

TODO: Shapes of first geometric order have a constant normal and do not need ξ.
"""
function getnormal(ξ, ::Type{T}, shapecoordinates, iface) where {T <: AbstractShape3D}
    n = parameterization_normal(ξ, T, shapecoordinates, iface)
    return n / norm(n)
end

"""
    parameterization_tangent(ξ, shapetype1D, shapecoordinates)
    parameterization_tangent(shapetype1D, jacobianNx1)

Return the tangent vector of the parameterization of a 1D shape in 2D or 3D space
at reference coordinates `ξ`.

NOTE: The returned vector does NOT have unit length so it can be used as the
differential tangent element in line integrals.

# References
https://en.wikipedia.org/wiki/Line_integral
"""
function parameterization_tangent(ξ, ::Type{T}, shapecoordinates) where {T <: AbstractShape1D}
    J = parameterization_jacobian(ξ, T, shapecoordinates)
    t = parameterization_tangent(T, J)
    return t
end

function parameterization_tangent(::Type{<:AbstractShape1D}, J::AbstractMatrix)
    @assert size(J, 2) == 1
    return J[:, 1]
end

"""
    parameterization_tangent(ξ, shapetype2D, shapecoordinates2D, iedge)
    parameterization_tangent(ξ, shapetype3D, shapecoordinates3D, iedge)
    parameterization_tangent(shapetype2D, jacobian2x2, iedge)
    parameterization_tangent(shapetype3D, jacobian3x3, iedge)

Return the tangent vector of edge `iedge` of the parameterization of a 2D or 3D shape
in 2D or 3D space at reference coordinates `ξ`. `ξ` must describe a point on the
local edge `iedge`.

WARNING: It is not checked if reference coordinates `ξ` actualle lie on edge `iedge`.

NOTE: The returned vector does NOT have unit length so it can be used as the
differential tangent element in line integrals.

# References
https://en.wikipedia.org/wiki/Line_integral
"""
function parameterization_tangent(ξ, T::Type{<:Union{AbstractShape2D, AbstractShape3D}}, shapecoordinates, iedge)
    J = parameterization_jacobian(ξ, T, shapecoordinates)
    return parameterization_tangent(T, J, iedge)
end

function parameterization_tangent(::Type{<:AbstractTriangle}, J::AbstractMatrix, iedge::Integer)
    @assert size(J) == (2, 2)
    iedge == 1 && return J[:, 1]
    iedge == 2 && return -J[:, 1] + J[:, 2]
    iedge == 3 && return -J[:, 2]
    throw(ArgumentError("no tangent defined for edge $iedge in AbstactTriangle"))
end

function parameterization_tangent(::Type{<:AbstractQuadrilateral}, J::AbstractMatrix, iedge::Integer)
    @assert size(J) == (2, 2)
    iedge == 1 && return  J[:, 1]
    iedge == 2 && return  J[:, 2]
    iedge == 3 && return -J[:, 1]
    iedge == 4 && return -J[:, 2]
    throw(ArgumentError("no tangent defined for edge $iedge in AbstractQuadrilateral"))
end

function parameterization_tangent(::Type{<:AbstractTetrahedron}, J::AbstractMatrix, iedge::Integer)
    @assert size(J) == (3, 3)
    iedge == 1 && return  J[:, 1]
    iedge == 2 && return -J[:, 1] + J[:, 2]
    iedge == 3 && return -J[:, 2]
    iedge == 4 && return  J[:, 3]
    iedge == 5 && return -J[:, 1] + J[:, 3]
    iedge == 6 && return -J[:, 2] + J[:, 2]
    throw(ArgumentError("no tangent defined for edge $iedge in AbstractTetrahedron"))
end

function parameterization_tangent(::Type{<:AbstractHexahedron}, J::AbstractMatrix, iedge::Integer)
    @assert size(J) == (3, 3)
    iedge ==  1 && return  J[:, 1]
    iedge ==  2 && return  J[:, 2]
    iedge ==  3 && return -J[:, 1]
    iedge ==  4 && return -J[:, 2]
    iedge ==  5 && return  J[:, 1]
    iedge ==  6 && return  J[:, 2]
    iedge ==  7 && return -J[:, 1]
    iedge ==  8 && return -J[:, 2]
    iedge ==  9 && return  J[:, 3]
    iedge == 10 && return  J[:, 3]
    iedge == 11 && return  J[:, 3]
    iedge == 12 && return  J[:, 3]
    throw(ArgumentError("no tangent defined for edge $iedge in AbstractHexahedron"))
end

function parameterization_tangent(::Type{<:AbstractWedge}, J::AbstractMatrix, iedge::Integer)
    @assert size(J) == (3, 3)
    iedge == 1 && return  J[:, 1]
    iedge == 2 && return -J[:, 1] + J[:, 2]
    iedge == 3 && return -J[:, 2]
    iedge == 4 && return  J[:, 1]
    iedge == 5 && return -J[:, 1] + J[:, 2]
    iedge == 6 && return -J[:, 2]
    iedge == 7 && return  J[:, 3]
    iedge == 8 && return  J[:, 3]
    iedge == 9 && return  J[:, 3]
    throw(ArgumentError("no tangent defined for edge $iedge in AbstractWedge"))
end

"""
    gettangent(ξ, shapetype1D, shapecoordinates3D)

Return the tangent vector (unit length) of the parameterization of a 1D shape
in 3D space at reference coordinates `ξ`.

TODO: Shapes of first geometric order have a constant normal and do not need ξ.
"""
function gettangent(ξ, ::Type{T}, shapecoordinates) where {T <: AbstractShape1D}
    t = parameterization_tangent(ξ, T, shapecoordinates)
    return t / norm(t)
end

"""
    gettangent(ξ, shapetype2D, shapecoordinates2D, iedge)
    gettangent(ξ, shapetype3D, shapecoordinates3D, iedge)

Return the tangent vector (unit length) of edge `iedge` of the parameterization
of a 2D shape in 2D space at reference coordinates `ξ`. `ξ` must describe a point on
the local edge `iedge`.

WARNING: It is not checked if reference coordinates `ξ` actualle lie on edge `iedge`.

TODO: Shapes of first geometric order have a constant tangent and do not need ξ.
"""
function gettangent(ξ, T::Type{<:Union{AbstractShape2D, AbstractShape3D}}, shapecoordinates, iedge)
    t = parameterization_tangent(ξ, T, shapecoordinates, iedge)
    return t / norm(t)
end

##############################
# Misc
##############################

struct Orientation
    flip::Bool
    ncircshift::Int
end

"""
    getorientation(shape::T) where T <: AbstractFace
    getorientation(shape::T) where T <: AbstractEdge

Return the `Orientation` of `shape` with respect to the global
orientation scheme.

The `Orientation` describes the number of
`flip` and `circshift` operations on `shape` in order to align it
with the global orientation specified by `to_global_orientation()`.

By convention, `flip()` must be performed before `circshift()`.

# See also:
`flip`
`circshift`
`to_global_orientation`
`flip_edgedofs_permutation`
`flip_facedofs_permutation`
`circshift_facedofs_permutation`

# Examples:
julia> shape = Quad8(2, 1, 3, 4, 5, 6, 7, 8);

julia> getorientation(shape)
Felder.Orientation(true, 2)
"""
getorientation(shape) = getorientation(shape, to_global_orientation(shape))

"""
    getorientation(a::T, b::T) where T <: AbstractFace
    getorientation(a::T, b::T) where T <: AbstractEdge

Return an `Orientation` containing (`flip`, `ncircshift`), where `flip` (`Bool`)
specifies whether a flip operation `flip(a)` is necessary and ncircshift is the
number of `circshift(a)` in order to align polygon shape `a` with `b`.

By convention, `flip()` must be performed before `circshift()`.

# Examples:
```julia-repl
julia> getorientation(Quad4(1,2,3,4), Quad4(1,2,3,4))
Orientation(false, 0)

julia> getorientation(Quad4(1,2,3,4), Quad4(3,4,1,2))
Orientation(false, 2)

julia> circshift(Quad4(1,2,3,4), 2) == Quad4(3,4,1,2)
true

julia> getorientation(Quad4(1,2,3,4), Quad4(3,2,1,4))
Orientation(true, 1)

julia> circshift(flip(Quad4(1,2,3,4)), 1) == Quad4(3,2,1,4)
true

julia> getorientation(Edge2(1,2), Edge2(1,2))
Orientation(false, 0)

julia> getorientation(Edge3(1,2,3), Edge3(2,1,3))
Orientation(true, 0)

julia> flip(Edge3(1,2,3)) == Edge3(2,1,3)
true
```
"""
function getorientation(a::T, b::T) where T <: AbstractFace
    doflip = false
    ncircshift = 0
    success = false

    a_ = deepcopy(a)

    for i in 1:nvertices(a_)
        if a_ == b
            success = true
            break
        end

        a_ = circshift(a_)
        ncircshift += 1
    end

    if !success
        a_ = flip(a)
        doflip = true
        ncircshift = 0

        for i in 1:nvertices(a_)
            if a_ == b
                success = true
                break
            end

            a_ = circshift(a_)
            ncircshift += 1
        end
    end

    success || throw(ArgumentError("Could not determine relative orientation of $a and $b."))

    return Orientation(doflip, ncircshift)
end

function getorientation(a::T, b::T) where T <: AbstractEdge
    if a == b
        return Orientation(false, 0)
    elseif flip(a) == b
        return Orientation(true, 0)
    else
        throw(ArgumentError("Could not determine relative orientation of $a and $b."))
    end
end

"""
    to_global_orientation(shape::AbstractEdge)
    to_global_orientation(shape::AbstractFace)

Return the `shape` with the nodes reordered according to a global orientation rule.
This is useful for identifying unique edges or faces in the mesh, or for vector basis
functions when a unique edge or face direction is required.

Note: The `shape` is returned with `shape.id` = 0.

# Examples:
```julia-repl
julia> to_global_orientation(Edge2(2,3))
Edge2([2, 3], 0)

julia> to_global_orientation(Edge2(3,2))
Edge2([2, 3], 0)

julia> to_global_orientation(Quad4(2,3,4,5))
Quad4([2, 3, 4, 5], 0)

julia> to_global_orientation(Quad4(4,3,2,5))
Quad4([2, 3, 4, 5], 0)
```
"""
function to_global_orientation(shape::Vertex)
    return Vertex(shape.n[1])
end

function to_global_orientation(shape::Edge2)
    n = shape.n
    if n[1] < n[2]
        return Edge2(n[1], n[2])
    else
        return Edge2(n[2], n[1])
    end
end

function to_global_orientation(shape::Edge3)
    n = shape.n
    if n[1] < n[2]
        return Edge3(n[1], n[2], n[3])
    else
        return Edge3(n[2], n[1], n[3])
    end
end

function to_global_orientation(shape::Tri3)
    n = shape.n
    imin = argmin((n[1], n[2], n[3]))
    if imin == 1
        if n[2] < n[3]
            return Tri3(n[1], n[2], n[3])
        else
            return Tri3(n[1], n[3], n[2])
        end
    elseif imin == 2
        if n[3] < n[1]
            return Tri3(n[2], n[3], n[1])
        else
            return Tri3(n[2], n[1], n[3])
        end
    else
        if n[1] < n[2]
            return Tri3(n[3], n[1], n[2])
        else
            return Tri3(n[3], n[2], n[1])
        end
    end
end

function to_global_orientation(shape::Tri6)
    n = shape.n
    imin = argmin((n[1], n[2], n[3]))
    if imin == 1
        if n[2] < n[3]
            return Tri6(n[1], n[2], n[3], n[4], n[5], n[6])
        else
            return Tri6(n[1], n[3], n[2], n[6], n[5], n[4])
        end
    elseif imin == 2
        if n[3] < n[1]
            return Tri6(n[2], n[3], n[1], n[5], n[6], n[4])
        else
            return Tri6(n[2], n[1], n[3], n[4], n[6], n[5])
        end
    else
        if n[1] < n[2]
            return Tri6(n[3], n[1], n[2], n[6], n[4], n[5])
        else
            return Tri6(n[3], n[2], n[1], n[5], n[4], n[6])
        end
    end
end

function to_global_orientation(shape::Quad4)
    n = shape.n
    imin = argmin((n[1], n[2], n[3], n[4]))
    if imin == 1
        if n[2] < n[4]
            return Quad4(n[1], n[2], n[3], n[4])
        else
            return Quad4(n[1], n[4], n[3], n[2])
        end
    elseif imin == 2
        if n[3] < n[1]
            return Quad4(n[2], n[3], n[4], n[1])
        else
            return Quad4(n[2], n[1], n[4], n[3])
        end
    elseif imin == 3
        if n[4] < n[2]
            return Quad4(n[3], n[4], n[1], n[2])
        else
            return Quad4(n[3], n[2], n[1], n[4])
        end
    else
        if n[1] < n[3]
            return Quad4(n[4], n[1], n[2], n[3])
        else
            return Quad4(n[4], n[3], n[2], n[1])
        end
    end
end

function to_global_orientation(shape::Quad8)
    n = shape.n
    imin = argmin((n[1], n[2], n[3], n[4]))
    if imin == 1
        if n[2] < n[4]
            return Quad8(n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8])
        else
            return Quad8(n[1], n[4], n[3], n[2], n[8], n[7], n[6], n[5])
        end
    elseif imin == 2
        if n[3] < n[1]
            return Quad8(n[2], n[3], n[4], n[1], n[6], n[7], n[8], n[5])
        else
            return Quad8(n[2], n[1], n[4], n[3], n[5], n[8], n[7], n[6])
        end
    elseif imin == 3
        if n[4] < n[2]
            return Quad8(n[3], n[4], n[1], n[2], n[7], n[8], n[5], n[6])
        else
            return Quad8(n[3], n[2], n[1], n[4], n[6], n[5], n[8], n[7])
        end
    else
        if n[1] < n[3]
            return Quad8(n[4], n[1], n[2], n[3], n[8], n[5], n[6], n[7])
        else
            return Quad8(n[4], n[3], n[2], n[1], n[7], n[6], n[5], n[8])
        end
    end
end

"""
    issameshape(shape1, shape2)

Return `true` if `shape1` and `shape2` are the same shape, i.e. the nodes describe
the same shape regardless of their specific order.

# Examples
```julia-repl
julia> issameshape(Vertex(2), Vertex(2))
true

julia> issameshape(Vertex(2), Vertex(3))
false

julia> issameshape(Tri3(2,3,4), Tri3(3,4,2))
true

julia> issameshape(Tri3(2,3,4), Tri3(3,4,6))
false

julia> issameshape(Tri3(2,3,4), Quad4(2,3,4,5))
false
```
"""
function issameshape(shape1::T, shape2::T) where {T <: AbstractShape0D}
    shape1.n == shape2.n
end

function issameshape(shape1::T, shape2::T) where {T <: AbstractShape1D}
    to_global_orientation(shape1) == to_global_orientation(shape2)
end

function issameshape(shape1::T, shape2::T) where {T <: AbstractShape2D}
    to_global_orientation(shape1) == to_global_orientation(shape2)
end

function issameshape(::T, ::T) where {T <: AbstractShape3D}
    error("not implemented for 3D shapes.")
end

issameshape(::AbstractShape, ::AbstractShape) = false
