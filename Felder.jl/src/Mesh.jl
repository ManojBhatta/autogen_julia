export AbstractFrame
export SpatialFrame
export MaterialFrame
export AbstractMesh, AbstractMesh1D, AbstractMesh2D, AbstractMesh3D
export Mesh1D, Mesh2D, Mesh3D
export Mesh
export elementtype, facettype
export ElementIndex, FacetIndex, EdgeIndex, VertexIndex

export ncoordinates, nvertices, nedges, nfaces, ncells
export nfacets, nelements
export getcoordinates
export getshapes
export getfacet, getelement
export getfacets, getelements
export getvertex, getedge, getface, getcell
export eachvertexgrouptag, eachedgegrouptag, eachboundarytag, eachdomaintag
export eachvertexgroup, eachedgegroup, eachboundary, eachdomain
export getvertexgroup, getedgegroup, getboundary, getdomain
export getneighbors
export getshape2color
export getcolor2elements
export eachelement
export eachelementindex
export eachfacetindex
export eachedgeindex
export eachvertexindex
export eachfacet
export eachedge
export eachvertex

#############################
# Frames
#############################

abstract type AbstractFrame end

"""
"""
struct SpatialFrame <: AbstractFrame end

"""
"""
struct MaterialFrame <: AbstractFrame end

##############################
# AbstractMesh
##############################

abstract type AbstractMesh{Dim} end

# Type aliases
const AbstractMesh1D = AbstractMesh{1}
const AbstractMesh2D = AbstractMesh{2}
const AbstractMesh3D = AbstractMesh{3}

Base.ndims(::AbstractMesh{Dim}) where {Dim} = Dim
Base.ndims(::Type{<:AbstractMesh{Dim}}) where {Dim} = Dim

ncoordinates(m::AbstractMesh) = length(m.coordinates)

nvertices(m::AbstractMesh1D) = length(m.facets)
nvertices(m::AbstractMesh2D) = length(m.vertices)
nvertices(m::AbstractMesh3D) = length(m.vertices)
nedges(m::AbstractMesh1D) = length(m.elements)
nedges(m::AbstractMesh2D) = length(m.facets)
nedges(m::AbstractMesh3D) = length(m.edges)
nfaces(m::AbstractMesh2D) = length(m.elements)
nfaces(m::AbstractMesh3D) = length(m.facets)
ncells(m::AbstractMesh3D) = length(m.elements)

nfacets(m::AbstractMesh) = length(m.facets)
nelements(m::AbstractMesh) = length(m.elements)

getshapes(m::AbstractMesh{N}, ::Val{N}) where {N} = m.elements
getshapes(m::AbstractMesh3D, ::Val{2}) = m.facets
getshapes(m::AbstractMesh3D, ::Val{1}) = m.edges
getshapes(m::AbstractMesh3D, ::Val{0}) = m.vertices
getshapes(m::AbstractMesh2D, ::Val{1}) = m.facets
getshapes(m::AbstractMesh2D, ::Val{0}) = m.vertices
getshapes(m::AbstractMesh1D, ::Val{0}) = m.facets

getfacet(m::AbstractMesh, i::Integer) = m.facets[i]
getelement(m::AbstractMesh, i::Integer) = m.elements[i]

getfacets(m::AbstractMesh) = m.facets
getelements(m::AbstractMesh) = m.elements

getvertex(m::AbstractMesh1D, i::Integer) = m.facets[i]
getvertex(m::AbstractMesh2D, i::Integer) = m.vertices[i]
getvertex(m::AbstractMesh3D, i::Integer) = m.vertices[i]
getedge(m::AbstractMesh1D, i::Integer) = m.elements[i]
getedge(m::AbstractMesh2D, i::Integer) = m.facets[i]
getedge(m::AbstractMesh3D, i::Integer) = m.edges[i]
getface(m::AbstractMesh2D, i::Integer) = m.elements[i]
getface(m::AbstractMesh3D, i::Integer) = m.facets[i]
getcell(m::AbstractMesh3D, i::Integer) = m.elements[i]

getcoordinates(m::AbstractMesh, ::SpatialFrame) = m.coordinates
getcoordinates(m::AbstractMesh, ::MaterialFrame) = m.initcoordinates

getvertexgroup(m::AbstractMesh, tag::Symbol)    = m.vertexgroups[tag]
getedgegroup(m::AbstractMesh, tag::Symbol)      = m.edgegroups[tag]
getboundary(m::AbstractMesh, tag::Symbol)       = m.boundaries[tag]
getdomain(m::AbstractMesh, tag::Symbol)         = m.domains[tag]

eachvertexgrouptag(m::AbstractMesh) = (tag for tag in keys(m.vertexgroups))
eachedgegrouptag(m::AbstractMesh) = (tag for tag in keys(m.edgegroups))
eachboundarytag(m::AbstractMesh) = (tag for tag in keys(m.boundaries))
eachdomaintag(m::AbstractMesh) = (tag for tag in keys(m.domains))

eachvertexgroup(m::AbstractMesh) = ((tag, group) for (tag, group) in m.vertexgroups)
eachedgegroup(m::AbstractMesh) = ((tag, group) for (tag, group) in m.edgegroups)
eachboundary(m::AbstractMesh) = ((tag, group) for (tag, group) in m.boundaries)
eachdomain(m::AbstractMesh) = ((tag, group) for (tag, group) in m.domains)

eachvertexgroup(m::AbstractMesh, tags) = ((tag, m.vertexgroups[tag]) for tag in tags)
eachedgegroup(m::AbstractMesh, tags) = ((tag, m.edgegroups[tag]) for tag in tags)
eachboundary(m::AbstractMesh, tags) = ((tag, m.boundaries[tag]) for tag in tags)
eachdomain(m::AbstractMesh, tags) = ((tag, m.domains[tag]) for tag in tags)

const MeshTagTypes = Union{Symbol, AbstractArray{Symbol}, NTuple{<:Any, Symbol}}

#############################
# Indices
#############################

abstract type AbstractMeshIndex end

"""
"""
struct VertexIndex <: AbstractMeshIndex
    shapeid::Int
    elementid::Int
    ilocal::Int
end

"""
"""
struct EdgeIndex <: AbstractMeshIndex
    shapeid::Int
    elementid::Int
    ilocal::Int
end

"""
"""
struct FacetIndex <: AbstractMeshIndex
    shapeid::Int
    elementid::Int
    ilocal::Int
end

"""
"""
struct ElementIndex <: AbstractMeshIndex
    shapeid::Int
end

@inline getelement(m::AbstractMesh, i::ElementIndex) = m.elements[i.shapeid]
@inline getelement(m::AbstractMesh, i::AbstractMeshIndex) = m.elements[i.elementid]

@inline Base.getindex(m::AbstractMesh, i::ElementIndex) = m.elements[i.shapeid]
@inline Base.getindex(m::AbstractMesh, i::FacetIndex) = m.facets[i.shapeid]
@inline Base.getindex(m::AbstractMesh3D, i::EdgeIndex) = m.edges[i.shapeid]
@inline Base.getindex(m::Union{AbstractMesh2D, AbstractMesh3D}, i::VertexIndex) = m.vertices[i.shapeid]

##############################
# Iterators
##############################

"""
"""
function eachelementindex(mesh::AbstractMesh, domaintags=eachdomaintag(mesh))
    (elementindex for tag in domaintags for elementindex in mesh.domains[tag])
end

eachelementindex(mesh::AbstractMesh, domain::Symbol) = eachelementindex(mesh, (domain,))

"""
"""
function eachelement(mesh::AbstractMesh, domaintags=eachdomaintag(mesh))
    (mesh[elementindex] for elementindex in eachelementindex(mesh, domaintags))
end

# ------------------------------------

"""
"""
function eachfacetindex(mesh::AbstractMesh, boundarytags)
    (facetindex for tag in boundarytags for facetindex in mesh.boundaries[tag])
end

eachfacetindex(mesh::AbstractMesh, boundarytag::Symbol) = eachfacetindex(mesh, (boundarytag,))

"""
"""
function eachfacet(mesh::AbstractMesh, boundarytags)
    (mesh[facetindex] for facetindex in eachfacetindex(mesh, boundarytags))
end

# ------------------------------------

"""
"""
function eachedgeindex(mesh::AbstractMesh3D, edgegrouptags)
    (edgeindex for tags in edgegrouptags for edgeindex in mesh.edgegroups[tags])
end

eachedgeindex(mesh::AbstractMesh3D, edgegrouptag::Symbol) = eachedgeindex(mesh, (edgegrouptag,))

"""
"""
function eachedge(mesh::AbstractMesh3D, edgegrouptags)
    (mesh[edgeindex] for edgeindex in eachedgeindex(mesh, edgegrouptags))
end

# ------------------------------------

"""
"""
function eachvertexindex(mesh::Union{AbstractMesh2D, AbstractMesh3D}, vertexgrouptags)
    (vertexindex for tag in vertexgrouptags for vertexindex in mesh.vertexgroups[tag])
end

eachvertexindex(mesh::Union{AbstractMesh2D, AbstractMesh3D}, vertexgrouptag::Symbol) = eachvertexindex(mesh, (vertexgrouptag,))

"""
"""
function eachvertex(mesh::Union{AbstractMesh2D, AbstractMesh3D}, vertexgrouptags)
    (mesh[vertexindex] for vertexindex in eachvertexindex(mesh, vertexgrouptags))
end

##############################
# Mesh
#############################

@with_kw struct Mesh1D{T1} <: AbstractMesh{1}
    coordinates::Vector{SVector{1, Float64}}
    elements::Vector{T1}
    facets::Vector{Vertex}                                  = Vector{Vertex}()

    domains::OrderedDict{Symbol, Vector{ElementIndex}}      = OrderedDict{Symbol,Vector{ElementIndex}}()
    boundaries::OrderedDict{Symbol, Vector{FacetIndex}}     = OrderedDict{Symbol,Vector{FacetIndex}}()

    domain2boundaries::OrderedDict{Symbol,Vector{Symbol}}   = OrderedDict{Symbol,Vector{Symbol}}()

    color2elements::Vector{Vector{ElementIndex}}            = Vector{ElementIndex}[]
    initcoordinates::Vector{SVector{1, Float64}}       = deepcopy(coordinates)
end

@with_kw struct Mesh2D{T1, T2} <: AbstractMesh{2}
    coordinates::Vector{SVector{2, Float64}}
    elements::Vector{T1}
    facets::Vector{T2}                                      = Vector{AbstractEdge}()
    vertices::Vector{Vertex}                                = Vector{Vertex}()

    domains::OrderedDict{Symbol, Vector{ElementIndex}}      = OrderedDict{Symbol,Vector{ElementIndex}}()
    boundaries::OrderedDict{Symbol, Vector{FacetIndex}}     = OrderedDict{Symbol,Vector{FacetIndex}}()
    vertexgroups::OrderedDict{Symbol, Vector{VertexIndex}}  = OrderedDict{Symbol,Vector{VertexIndex}}()

    domain2boundaries::OrderedDict{Symbol,Vector{Symbol}}   = OrderedDict{Symbol,Vector{Symbol}}()
    domain2vertexgroups::OrderedDict{Symbol,Vector{Symbol}} = OrderedDict{Symbol,Vector{Symbol}}()

    color2elements::Vector{Vector{ElementIndex}}            = Vector{ElementIndex}[]
    initcoordinates::Vector{SVector{2, Float64}}       = deepcopy(coordinates)
end

@with_kw struct Mesh3D{T1, T2, T3} <: AbstractMesh{3}
    coordinates::Vector{SVector{3, Float64}}
    elements::Vector{T1}
    facets::Vector{T2}                                      = Vector{AbstractFace}()
    edges::Vector{T3}                                       = Vector{AbstractEdge}()
    vertices::Vector{Vertex}                                = Vector{Vertex}()

    domains::OrderedDict{Symbol, Vector{ElementIndex}}      = OrderedDict{Symbol,Vector{ElementIndex}}()
    boundaries::OrderedDict{Symbol, Vector{FacetIndex}}     = OrderedDict{Symbol,Vector{FacetIndex}}()
    edgegroups::OrderedDict{Symbol, Vector{EdgeIndex}}      = OrderedDict{Symbol,Vector{EdgeIndex}}()
    vertexgroups::OrderedDict{Symbol, Vector{VertexIndex}}  = OrderedDict{Symbol,Vector{VertexIndex}}()

    domain2boundaries::OrderedDict{Symbol,Vector{Symbol}}   = OrderedDict{Symbol,Vector{Symbol}}()
    domain2edgegroups::OrderedDict{Symbol,Vector{Symbol}}   = OrderedDict{Symbol,Vector{Symbol}}()
    domain2vertexgroups::OrderedDict{Symbol,Vector{Symbol}} = OrderedDict{Symbol,Vector{Symbol}}()

    color2elements::Vector{Vector{ElementIndex}}            = Vector{ElementIndex}[]
    initcoordinates::Vector{SVector{3, Float64}}       = deepcopy(coordinates)
end

Base.broadcastable(m::AbstractMesh) = Ref(m)

elementtype(::Mesh1D{T1}) where {T1} = T1
elementtype(::Mesh2D{T1}) where {T1} = T1
elementtype(::Mesh3D{T1}) where {T1} = T1

facettype(::Mesh1D) = Vertex
facettype(::Mesh2D{T1, T2}) where {T1, T2} = T2
facettype(::Mesh3D{T1, T2}) where {T1, T2} = T2

function Base.show(io::Core.IO, ::MIME"text/plain", m::AbstractMesh1D)
    println(io, typeof(m), " with:")
    println(io, "  Element type (edges):    ", eltype(getelements(m)))
    println(io, "  Facet type (vertices):   ", eltype(getfacets(m)))
    println(io, "  Number of coordinates:   ", ncoordinates(m))
    println(io, "  Number of elements:      ", nelements(m))
    println(io, "  Number of facets:        ", nfacets(m))
    println(io, "  Number of colors:        ", length(m.color2elements))
    println(io, "  Size in memory:          ", Base.format_bytes(Base.summarysize(m)))
    println(io, "  Domains:                 ", collect(eachdomaintag(m)))
    print(  io, "  Boundaries:              ", collect(eachboundarytag(m)))
end

function Base.show(io::Core.IO, ::MIME"text/plain", m::AbstractMesh2D)
    println(io, typeof(m), " with:")
    println(io, "  Element type (faces):    ", eltype(getelements(m)))
    println(io, "  Facet type (edges):      ", eltype(getfacets(m)))
    println(io, "  Number of coordinates:   ", ncoordinates(m))
    println(io, "  Number of elements:      ", nelements(m))
    println(io, "  Number of facets:        ", nfacets(m))
    println(io, "  Number of vertices:      ", nvertices(m))
    println(io, "  Number of colors:        ", length(m.color2elements))
    println(io, "  Size in memory:          ", Base.format_bytes(Base.summarysize(m)))
    println(io, "  Domains:                 ", collect(eachdomaintag(m)))
    println(io, "  Boundaries:              ", collect(eachboundarytag(m)))
    print(  io, "  Vertex groups:           ", collect(eachvertexgrouptag(m)))
end

function Base.show(io::Core.IO, ::MIME"text/plain", m::AbstractMesh3D)
    println(io, typeof(m), " with:")
    println(io, "  Element type (cells):    ", eltype(getelements(m)))
    println(io, "  Facet type (faces):      ", eltype(getfacets(m)))
    println(io, "  Edge type (edges):       ", eltype(m.edges))
    println(io, "  Number of coordinates:   ", ncoordinates(m))
    println(io, "  Number of elements:      ", nelements(m))
    println(io, "  Number of facets:        ", nfacets(m))
    println(io, "  Number of edges:         ", nedges(m))
    println(io, "  Number of vertices:      ", nvertices(m))
    println(io, "  Number of colors:        ", length(m.color2elements))
    println(io, "  Size in memory:          ", Base.format_bytes(Base.summarysize(m)))
    println(io, "  Domains:                 ", collect(eachdomaintag(m)))
    println(io, "  Boundaries:              ", collect(eachboundarytag(m)))
    println(io, "  Edge groups:             ", collect(eachedgegrouptag(m)))
    print(  io, "  Vertex groups:           ", collect(eachvertexgrouptag(m)))
end

##############################
# Mesh Constructors
#############################

"""
    m = Mesh(filename)
    m = Mesh1D(filename)
    m = Mesh2D(filename)
    m = Mesh3D(filename)

Construct mesh from file specified by `filename`.

The constructor `Mesh()` detects how many dimensions the mesh points
span and returns an instance of `Mesh1D`. `Mesh2D` or `Mesh3D`.

The mesh can be forced to be constructed in higher dimensions than the
span of the mesh points (e.g. 2D plane in 3D space) by explicitly using
the constructors `Mesh1D()`, `Mesh2D()` or `Mesh3D()`.

# Example
```julia-repl
julia> m = Mesh("test/meshes/mesh1D.unv"); typeof(m)
Mesh1D{Edge2{Int64}}

julia> m = Mesh("test/meshes/mesh2D.unv"); typeof(m)
Mesh2D{Tri3{Int64}, Edge2{Int64}}

julia> m = Mesh("test/meshes/mesh3D.unv"); typeof(m)
Mesh3D{Tet4{Int64}, Tri3{Int64}, Edge2{Int64}}

julia> m = Mesh3D("test/meshes/mesh2D.unv"); typeof(m)
Mesh3D{AbstractVolume, Tri3{Int64}, Edge2{Int64}}
```
"""
function Mesh(filename::String; kwargs...)
    fileext = lowercase(splitext(filename)[end])
    if fileext == ".unv"
        mesh = Mesh(UnvMesh(filename); kwargs...)
    else
        error("Import of mesh file type $fileext not implemented!")
    end
    checkgrouptags(mesh)
    return mesh
end

##############################
# Misc
##############################

"""
    getboundingbox(mesh1D, frame=SpatialFrame()) -> (xmin, xmax)
    getboundingbox(mesh2D, frame=SpatialFrame()) -> (xmin, xmax), (ymin, ymax)
    getboundingbox(mesh3D, frame=SpatialFrame()) -> (xmin, xmax), (ymin, ymax), (zmin, zmax)

Return the bounding box of the mesh in the specified frame (default: `SpatialFrame`).

# Example
```julia-repl
julia> getboundingbox(TestMesh1D())
(0.0, 2.0)

julia> getboundingbox(TestMesh2D())
((0.0, 2.0), (0.0, 1.0))

julia> getboundingbox(TestMesh3D(), MaterialFrame())
((0.0, 2.0), (0.0, 1.0), (0.0, 0.5))
```
"""
getboundingbox(m::AbstractMesh1D, frame=SpatialFrame()) = getboundingbox1D(getcoordinates(m, frame))
getboundingbox(m::AbstractMesh2D, frame=SpatialFrame()) = getboundingbox2D(getcoordinates(m, frame))
getboundingbox(m::AbstractMesh3D, frame=SpatialFrame()) = getboundingbox3D(getcoordinates(m, frame))

"""
    getnodeconnectivity(shapes)

Return the node indices each node is "connected" to (including itself).
"""
function getnodeconnectivity(shapes::Vector{<:AbstractShape})
    shape2nodes = map(getnodes, shapes)
    node2shapes = invert_map(shape2nodes)

    nodeconnectivity = [Int[] for i in eachindex(node2shapes)]
    for i in eachindex(node2shapes)
        for e in node2shapes[i]
            for j in shape2nodes[e]
                unique_sorted_insert!(nodeconnectivity[i], j)
            end
        end
    end

    return nodeconnectivity
end

"""
    getneighbors(shapes) -> Vector{Vector{Int}}

Return the indices of neighboring shapes for all `shapes`.
Neighboring shapes share at least one vertex. (including itself)

Use `getneighbors(shapes, i)` if only the neighbors of
a single shape is needed.
"""
function getneighbors(shapes::Vector{<:AbstractShape})
    shape2vertices = map(getvertexnodes, shapes)
    vertex2shapes = invert_map(shape2vertices)

    shapeneighbors = [Int[] for i in eachindex(shape2vertices)]
    for i in eachindex(shape2vertices)
        for vertex in shape2vertices[i]
            for j in vertex2shapes[vertex]
                unique_sorted_insert!(shapeneighbors[i], j)
            end
        end
    end

    return shapeneighbors
end

"""
    getneighbors(shapes, i) -> Vector{Int}

Return the indices of neighboring shapes for `shapes[i]`.
Neighboring shapes share at least one vertex. (including itself)

Use `getneighbors(shapes)` to get the neighbors of all shapes.
"""
function getneighbors(shapes::Vector{<:AbstractShape}, i::Integer)
    shape2vertices = map(getvertexnodes, shapes)
    vertex2shapes = invert_map(shape2vertices)

    shapeneighbors = Int[]
    for vertex in shape2vertices[i]
        for j in vertex2shapes[vertex]
            unique_sorted_insert!(shapeneighbors, j)
        end
    end

    return shapeneighbors
end

"""
    find_a_parent(shape, elements)

Return the first parent of a vertex, edge or face as a integer tuple `(e, i)`,
where `e` is the element index and `i` the local index to the vertex, edge or face
in that element.

Only the first parent is returned which is found by iterating over all elements
connected to the `shape` vertices.

If the node2elements mapping is already known, it can be passed as an
optional argument to avoid redundant calculations.
"""
function find_a_parent(vertex::AbstractShape0D, elements::Vector{<:AbstractShape},
                       node2elements::Vector{<:Vector{<:Integer}}=invert_map(getvertexnodes, elements))
    for e in node2elements[first(getnodes(vertex))]
        for (i, local_vertex) in enumerate(getlocalvertices(elements[e]))
            if issameshape(vertex, local_vertex)
                return e, i
            end
        end
    end
    throw(ArgumentError("could not find a local $vertex in elements."))
end

function find_a_parent(edge::AbstractShape1D, elements::Vector{<:AbstractShape},
                       node2elements::Vector{<:Vector{<:Integer}}=invert_map(getvertexnodes, elements))
    for e in node2elements[first(getvertexnodes(edge))]
        for (i, local_edge) in enumerate(getlocaledges(elements[e]))
            if issameshape(edge, local_edge)
                return e, i
            end
        end
    end
    throw(ArgumentError("could not find a local $edge in elements."))
end

function find_a_parent(face::AbstractShape2D, elements::Vector{<:AbstractShape},
                       node2elements::Vector{<:Vector{<:Integer}}=invert_map(getvertexnodes, elements))
    for e in node2elements[first(getvertexnodes(face))]
        for (i, local_face) in enumerate(getlocalfaces(elements[e]))
            if issameshape(face, local_face)
                return e, i
            end
        end
    end
    throw(ArgumentError("could not find a local $face in elements."))
end

"""
Return a generator `(e, i)` over all parents.
"""
function iterate_all_parents(vertex::AbstractShape0D, elements::Vector{<:AbstractShape},
                       node2elements::Vector{<:Vector{<:Integer}}=invert_map(getvertexnodes, elements))
    ((e, i) for e in node2elements[first(getnodes(vertex))]
                for (i, local_vertex) in enumerate(getlocalvertices(elements[e]))
                    if issameshape(vertex, local_vertex))
end

function iterate_all_parents(edge::AbstractShape1D, elements::Vector{<:AbstractShape},
                       node2elements::Vector{<:Vector{<:Integer}}=invert_map(getvertexnodes, elements))
    ((e, i) for e in node2elements[first(getvertexnodes(edge))]
                for (i, local_edge) in enumerate(getlocaledges(elements[e]))
                    if issameshape(edge, local_edge))
end

function iterate_all_parents(face::AbstractShape2D, elements::Vector{<:AbstractShape},
                       node2elements::Vector{<:Vector{<:Integer}}=invert_map(getvertexnodes, elements))
    ((e, i) for e in node2elements[first(getvertexnodes(face))]
                for (i, local_face) in enumerate(getlocalfaces(elements[e]))
                    if issameshape(face, local_face))
end

function assign_colors_greedy(elementneighbors::Vector{<:AbstractVector{<:Integer}})
    element2color = zeros(Int, length(elementneighbors))

    max_color = 0
    for i in eachindex(elementneighbors)
        neighborcolors = @view element2color[elementneighbors[i]]
        for c in 1:(max_color + 1)
            if c ∉ neighborcolors
                element2color[i] = c
                max_color = max(max_color, c)
                break
            end
        end
    end

    return element2color
end

function getshape2color(shapes::Vector{<:AbstractShape})
    shapeneighbors = getneighbors(shapes)
    shape2color = assign_colors_greedy(shapeneighbors)
    return shape2color
end

function getcolor2elements(shapes::Vector{<:AbstractShape})
    element2color = getshape2color(shapes)
    color2elements = getcolor2elements(element2color)
    return color2elements
end

function getcolor2elements(element2color::Vector{<:Integer})
    max_color = maximum(element2color)
    color2elements = [ElementIndex[] for i in 1:max_color]

    for (i, c) in enumerate(element2color)
        push!(color2elements[c], ElementIndex(i))
    end

    return color2elements
end

"""
    renamegroups!(mesh, oldname => newname)
    renamegroups!(mesh, oldname => newname, oldname => newname)
    renamegroups!(mesh, Dict(oldname => newname, ...))

Rename a group in `mesh`.
"""
function renamegroups!(mesh::AbstractMesh, names::AbstractDict{Symbol, Symbol})
    if mesh isa AbstractMesh1D
        groupfields = (mesh.domains, mesh.boundaries)
        subgroupfields = (mesh.domain2boundaries,)
    elseif mesh isa AbstractMesh2D
        groupfields = (mesh.domains, mesh.boundaries, mesh.vertexgroups)
        subgroupfields = (mesh.domain2boundaries, mesh.domain2vertexgroups)
    else
        groupfields = (mesh.domains, mesh.boundaries, mesh.edgegroups, mesh.vertexgroups)
        subgroupfields = (mesh.domain2boundaries, mesh.domain2edgegroups, mesh.domain2vertexgroups)
    end

    for (oldname, newname) in names
        found = false
        for groups in groupfields
            if oldname in keys(groups)
                replacekey!(groups, oldname, newname)
                found = true
            end
        end
        if !found
            throw(KeyError(oldname))
        end
    end

    for (oldname, newname) in names
        for subgroups in subgroupfields
            for subgroupnames in values(subgroups)
                replace!(subgroupnames, oldname => newname)
            end
            if oldname in keys(subgroups)
                replacekey!(subgroups, oldname, newname)
            end
        end
    end
    return mesh
end

renamegroups!(mesh::AbstractMesh, names::Pair...) = renamegroups!(mesh, OrderedDict(names))
renamegroups!(mesh::AbstractMesh, names::NTuple{<:Any, <:Pair}) = renamegroups!(mesh, OrderedDict(names))

function checkgrouptags(mesh::AbstractMesh)
    tags = Symbol[]
    append!(tags, keys(mesh.domains))
    append!(tags, keys(mesh.boundaries))
    hasproperty(mesh, :edgegroups)   && append!(tags, keys(mesh.edgegroups))
    hasproperty(mesh, :vertexgroups) && append!(tags, keys(mesh.vertexgroups))
    sort!(tags)
    for i in 2:length(tags)
        if tags[i] == tags[i-1]
            throw(ArgumentError("duplicate group tag: $tags[i]. Group tags must be unique."))
        end
    end
end

function findsubgroups(domains, boundaries)
    d = OrderedDict{Symbol, Vector{Symbol}}()
    for (elementgroupkey, elementgroup) in domains
        syms = Symbol[]
        for (facetgroupkey, facetgroup) in boundaries
            if containsgroup(elementgroup, facetgroup)
                push!(syms, facetgroupkey)
            end
        end
        d[elementgroupkey] = syms
    end
    return d
end

function containsgroup(elementindices, facegroup)
    nelements = maximum(x -> x.shapeid, elementindices)
    x = falses(nelements)
    for elementindex in elementindices
        x[elementindex.shapeid] = true
    end
    for faceindex in facegroup
        faceindex.elementid <= nelements || continue
        if x[faceindex.elementid]
            return true
        end
    end
    return false
end
