export DofMap
export MeshDofs
export getconnectivity
export renumber_dofs!
export map_nodes2dofs

##########################
# Dof Assignment
##########################

assigndofs(shapes, args...; kwargs...) = assigndofs(collect(shapes), args...; kwargs...)

#--------------------------
# Lagrange Shape Functions
#--------------------------

function assigndofs(shapes::Vector{<:AbstractShape1D}, shapefunc::AbstractLagrange; renumber=true)
    node2shapes = invert_map(getvertexnodes, shapes)

    shape2vertices, nvertices = assign_internal_vertices(shapes, node2shapes)

    # Global orientation
    vertexdofs = [Int[] for _ in 1:nvertices]
    edgedofs   = [Int[] for _ in 1:length(shapes)]
    counter = 0
    for (i, shape) in enumerate(shapes)
        shapetype = typeof(shape)
        for (j, n) in zip(shape2vertices[i], ndofs_vertex(shapetype, shapefunc))
            _, counter = _assign_n_dofs!(vertexdofs[j], n, counter)
        end
        _, counter = _assign_n_dofs!(edgedofs[i], ndofs_edge(shapetype, shapefunc)[1], counter)
    end

    # Append in order: vertices -> edge
    shapedofs = [Int[] for _ in 1:length(shapes)]
    for (i, shape) in enumerate(shapes)
        for v in shape2vertices[i]
            append!(shapedofs[i], vertexdofs[v])
        end
        append!(shapedofs[i], edgedofs[i])
    end

    if renumber
        renumber_dofs!(shapedofs)
    end

    return shapedofs, counter
end

function assigndofs(shapes::Vector{<:AbstractShape2D}, shapefunc::AbstractLagrange; renumber=true)
    node2shapes = invert_map(getvertexnodes, shapes)

    shape2vertices, nvertices = assign_internal_vertices(shapes, node2shapes)
    shape2edges,    nedges    = assign_internal_edges(shapes, node2shapes)

    # Global orientation
    vertexdofs = [Int[] for _ in 1:nvertices]
    edgedofs   = [Int[] for _ in 1:nedges]
    facedofs   = [Int[] for _ in 1:length(shapes)]
    counter = 0
    for (i, shape) in enumerate(shapes)
        shapetype = typeof(shape)
        for (j, n) in zip(shape2vertices[i], ndofs_vertex(shapetype, shapefunc))
            _, counter = _assign_n_dofs!(vertexdofs[j], n, counter)
        end
        for (j, n) in zip(shape2edges[i], ndofs_edge(shapetype, shapefunc))
            _, counter = _assign_n_dofs!(edgedofs[j], n, counter)
        end
        _, counter = _assign_n_dofs!(facedofs[i], ndofs_face(shapetype, shapefunc)[1], counter)
    end

    # Append in order: vertices -> edges -> face
    shapedofs = [Int[] for _ in 1:length(shapes)]
    for (i, shape) in enumerate(shapes)
        for v in shape2vertices[i]
            append!(shapedofs[i], vertexdofs[v])
        end
        for (j, e) in enumerate(shape2edges[i])
            localedge = getlocaledge(shape, j)
            append_dofs_in_local_orientation!(shapedofs[i], edgedofs[e], localedge, shapefunc)
        end
        append!(shapedofs[i], facedofs[i])
    end

    if renumber
        renumber_dofs!(shapedofs)
    end

    return shapedofs, counter
end

function assigndofs(shapes::Vector{<:AbstractShape3D}, shapefunc::AbstractLagrange; renumber=true)
    node2shapes = invert_map(getvertexnodes, shapes)

    shape2vertices, nvertices = assign_internal_vertices(shapes, node2shapes)
    shape2edges,    nedges    = assign_internal_edges(shapes, node2shapes)
    shape2faces,    nfaces    = assign_internal_faces(shapes, node2shapes)

    # Global orientation
    vertexdofs = [Int[] for _ in 1:nvertices]
    edgedofs   = [Int[] for _ in 1:nedges]
    facedofs   = [Int[] for _ in 1:nfaces]
    celldofs   = [Int[] for _ in 1:length(shapes)]
    counter = 0
    for (i, shape) in enumerate(shapes)
        shapetype = typeof(shape)
        for (j, n) in zip(shape2vertices[i], ndofs_vertex(shapetype, shapefunc))
            _, counter = _assign_n_dofs!(vertexdofs[j], n, counter)
        end
        for (j, n) in zip(shape2edges[i], ndofs_edge(shapetype, shapefunc))
            _, counter = _assign_n_dofs!(edgedofs[j], n, counter)
        end
        for (j, n) in zip(shape2faces[i], ndofs_face(shapetype, shapefunc))
            _, counter = _assign_n_dofs!(facedofs[j], n, counter)
        end
        _, counter = _assign_n_dofs!(celldofs[i], ndofs_cell(shapetype, shapefunc)[1], counter)
    end

    # Append in order: vertices -> edges -> faces -> cell
    shapedofs = [Int[] for _ in 1:length(shapes)]
    for (i, shape) in enumerate(shapes)
        for v in shape2vertices[i]
            append!(shapedofs[i], vertexdofs[v])
        end
        for (j, e) in enumerate(shape2edges[i])
            localedge = getlocaledge(shape, j)
            append_dofs_in_local_orientation!(shapedofs[i], edgedofs[e], localedge, shapefunc)
        end
        for (j, f) in enumerate(shape2faces[i])
            localface = getlocalface(shape, j)
            append_dofs_in_local_orientation!(shapedofs[i], facedofs[f], localface, shapefunc)
        end
        append!(shapedofs[i], celldofs[i])
    end

    if renumber
        renumber_dofs!(shapedofs)
    end

    return shapedofs, counter
end

function _assign_n_dofs!(dofs::Vector{Int}, n, counter)
    if isempty(dofs)
        for _ in 1:n
            counter += 1
            push!(dofs, counter)
        end
    else
        if length(dofs) ≠ n
            error("tried to assign $n dofs to index $i, but $(length(dofs)) dofs are already assigned")
        end
    end
    return dofs, counter
end

function append_dofs_in_local_orientation!(collection, globaldofs, shape, shapefunc)
    perm_ind = local_dof_permutation(shape, shapefunc)

    @assert length(perm_ind) == length(globaldofs)
    for i in perm_ind
        push!(collection, globaldofs[i])
    end

    return collection
end

function local_dof_permutation(shape::AbstractShape1D, shapefunc)
    orientation = getorientation(shape)

    flip_perm = flip_edgedofs_permutation(typeof(shape), shapefunc)

    if orientation.flip
        return flip_perm
    else
        return typeof(flip_perm)(1:length(flip_perm))
    end
end

function local_dof_permutation(shape::AbstractShape2D, shapefunc)
    orientation = getorientation(shape)

    flip_perm = flip_facedofs_permutation(typeof(shape), shapefunc)
    circshift_perm = circshift_facedofs_permutation(typeof(shape), shapefunc)

    perm_ind = typeof(flip_perm)(1:length(flip_perm))

    if orientation.flip
        perm_ind = perm_ind[flip_perm]
    end

    for _ in 1:orientation.ncircshift
        perm_ind = perm_ind[circshift_perm]
    end

    return perm_ind
end

#--------------------------
# Constant Shape Functions
#--------------------------

function assigndofs(shapes::Vector{<:AbstractShape}, ::ConstantShapeFunction; renumber=true)
    shapedofs = Vector{Int}[]
    resize!(shapedofs, length(shapes))

    counter = 0
    for i in eachindex(shapes)
        counter += 1
        shapedofs[i] = [counter]
    end

    if renumber
        # Ignored
    end

    return shapedofs, counter
end

###############################################
# Assign internal faces/edges/vertices indices
###############################################

function assign_internal_vertices(shapes::Vector{<:AbstractShape}, node2shapes=invert_map(getvertexnodes, shapes))
    shape2vertices = [zeros(Int, nvertices(typeof(shape))) for shape in shapes]
    counter = 0
    for (i_shape, shape) in enumerate(shapes)
        for i_vertex in 1:nvertices(typeof(shape))
            # Check if this face already been numbered. If so, skip.
            shape2vertices[i_shape][i_vertex] == 0 || continue

            # Create new face
            counter += 1
            shape2vertices[i_shape][i_vertex] = counter

            # Assign edge index to connected elements (i.e. mark as visited)
            vertex = getlocalvertex(shape, i_vertex)
            for (e, i) in iterate_all_parents(vertex, shapes, node2shapes)
                shape2vertices[e][i] = counter
            end
        end
    end
    return shape2vertices, counter
end

function assign_internal_edges(shapes::Vector{<:AbstractShape}, node2shapes=invert_map(getvertexnodes, shapes))
    shape2edges = [zeros(Int, nedges(typeof(shape))) for shape in shapes]
    counter = 0
    for (i_shape, shape) in enumerate(shapes)
        for i_edge in 1:nedges(typeof(shape))
            # Check if this face already been numbered. If so, skip.
            shape2edges[i_shape][i_edge] == 0 || continue

            # Create new face
            counter += 1
            shape2edges[i_shape][i_edge] = counter

            # Assign edge index to connected elements (i.e. mark as visited)
            localedge = getlocaledge(shape, i_edge)
            for (e, i) in iterate_all_parents(localedge, shapes, node2shapes)
                shape2edges[e][i] = counter
            end
        end
    end
    return shape2edges, counter
end

nfaces(::Type{<:AbstractShape1D}) = 0

function assign_internal_faces(shapes::Vector{<:AbstractShape}, node2shapes=invert_map(getvertexnodes, shapes))
    shape2faces = [zeros(Int, nfaces(typeof(shape))) for shape in shapes]
    counter = 0
    for (i_shape, shape) in enumerate(shapes)
        for i_face in 1:nfaces(typeof(shape))
            # Check if this face already been numbered. If so, skip.
            shape2faces[i_shape][i_face] == 0 || continue

            # Create new face
            counter += 1
            shape2faces[i_shape][i_face] = counter

            # Assign face index to connected elements (i.e. mark as visited)
            localface = getlocalface(shape, i_face)
            for (e, i) in iterate_all_parents(localface, shapes, node2shapes)
                shape2faces[e][i] = counter
            end
        end
    end
    return shape2faces, counter
end

##########################
# DofMap
##########################

abstract type AbstractDofMap end

struct DofMap{Sf} <: AbstractDofMap
    shapefunc::Sf
    tags::Vector{Symbol}
    ndofs::Int
    dofs::Vector{Vector{Int}}
end

sftype(::DofMap{Sf}) where {Sf} = Sf
ndofs(dofmap::AbstractDofMap) = dofmap.ndofs
@inline Base.getindex(dofmap::AbstractDofMap, i::Integer) = dofmap.dofs[i]

"""
    DofMap()

Constructs a degree-of-freedom map (DofMap) for a given shape function and mesh. The DofMap assigns a unique
index to each degree of freedom in the mesh, which is used to assemble the global system of equations.

# Examples
"""
function DofMap(shapefunc::AbstractShapeFunctions, shapes; renumber=true)
    dofs, n_dofs = assigndofs(shapes, shapefunc; renumber)
    return DofMap{typeof(shapefunc)}(shapefunc, Symbol[], n_dofs, dofs)
end

##########################
# Node renumbering
##########################

"""
    renumber_dofs!(dofmap)
    renumber_dofs!(shapedofs::Vector{Vector{Int}})

Renumber the degrees of freedom using the reverse Cuthill-McKee algorithm
to reduce the bandwidth of the assembled system matrix.

# See also
cuthill_mckee_renumbering

# Example
```julia-repl
julia> dofmap = DofMap(Mesh("test/meshes/mesh3D_Tet10.unv"), Lagrange(3), renumber=false);

julia> preassemble(dofmap.dofs)
2223×2223 SparseMatrixCSC{Float64, Int64} with 90355 stored entries:
⣿⣿⣿⣿⣿⣿⣿⣿⡧⣭⠭⡭⠅⣭⣭⣭
⣿⣿⣿⣿⣿⣿⣿⣿⡷⠯⠱⠂⠆⠺⠯⠿
⣿⣿⣿⣿⣿⣿⣿⣿⡿⣏⣍⣭⡇⡽⣋⡛
⣿⣿⣿⣿⣿⣿⣿⣿⣗⣿⡙⡒⡃⢛⣟⣟
⡍⣯⡽⡏⡿⢯⣽⣽⣿⣿⣿⣿⣿⣿⣿⣿
⡇⡧⠱⠂⡇⣽⢳⠨⣿⣿⣿⣿⣿⣿⣿⣿
⡅⣥⣨⡁⣍⡭⣭⢈⣿⣿⣿⣿⣿⣿⣿⣿
⡇⣿⣯⡇⣯⠸⣿⢽⣿⣿⣿⣿⣿⣿⣿⣿

julia> renumber_dofs!(dofmap);

julia> preassemble(dofmap.dofs)
2223×2223 SparseMatrixCSC{Float64, Int64} with 90355 stored entries:
⣿⣿⣿⣦⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠻⣿⣿⣿⣿⣦⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠘⠻⣿⢿⣷⣶⣿⢦⣄⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠹⣼⣿⡿⣯⣿⣿⣷⣄⠀⠀⠀⠀
⠀⠀⠀⠀⠈⢷⣿⣿⣿⣿⣿⢻⣳⡄⠀⠀
⠀⠀⠀⠀⠀⠀⠙⢿⣿⣛⣻⣾⣿⣓⣦⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⠾⢿⢻⣿⣿⣿⡇
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⠿⠿⣿⣿
```
"""
function renumber_dofs!(shapedofs::Vector{<:AbstractVector{Int}})
    connectivity = getconnectivity(shapedofs)
    renumbermap = cuthill_mckee_renumbering(connectivity)
    for dofs in shapedofs
        for i in eachindex(dofs)
            dofs[i] = renumbermap[dofs[i]]
        end
    end
    return shapedofs
end

renumber_dofs!(dofmap::AbstractDofMap) = renumber_dofs!(dofmap.dofs)

"""
    cuthill_mckee_renumbering(connectivity)

Node numbering algorithm based on the reverse Cuthill–McKee algorithm. Vector
`connectivity[i]` contains the indices of nodes connected to node `i`.

https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
http://www.juliafem.org/examples/2017-08-29-reordering-nodes-with-the-RCM-algorithm
"""
function cuthill_mckee_renumbering(connectivity::Vector{<:AbstractVector{Int}})
    visited = falses(length(connectivity))

    # Start from minimum degree node (ignore empty nodes)
    n, index = findmin(x -> isempty(x) ? Inf : length(x), connectivity)

    R = [index]
    visited[index] = true

    A_i = Int[]
    A_length = Int[]
    perm = Int[]

    i = 0
    while i < length(R)
        i += 1
        empty!(A_i)
        empty!(A_length)

        for j in connectivity[R[i]]
            if !(visited[j])
                push!(A_i, j)
                push!(A_length, length(connectivity[j]))
            end
        end

        length(A_i) > 0 || continue

        resize!(perm, length(A_length))
        sortperm!(perm, A_length)

        for j in perm
            push!(R, A_i[j])
            visited[A_i[j]] = true
        end
    end

    reverse!(R) # REVERSE Cuthill–McKee
    renumbermap = invert_injective_map(R)

    return renumbermap
end

##############################
# Misc
##############################

"""
    getconnectivity(shapedofs::Vector{<:AbstractVector{Int}})

Return the dofs each dof is connected to via adjacent elements, including
itself.
"""
function getconnectivity(shapedofs)
    dof2element = invert_map(shapedofs)

    connectivity = [Int[] for _ in eachindex(dof2element)]
    for i in eachindex(connectivity)
        for e in dof2element[i]
            for j in shapedofs[e]
                unique_sorted_insert!(connectivity[i], j)
            end
        end
    end

    return connectivity
end
