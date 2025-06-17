export TestMesh1D
export TestMesh2D
export TestMesh3D
export linsteprange
export AnnularSectorMesh
export AnnularSectorMesh2D
export AnnularSectorMesh3D
export SphericalShellPanelMesh

"""
Mesh group names on a 1D structured grid.

    ○─────> x

2 boundaries:

    ○─────────────○
    b1            b2

"""
meshgroupnames1D = (
    :domain_1   => :domain_1,
    :boundary_1 => :boundary_left,
    :boundary_2 => :boundary_right,
)

"""
Mesh group names on a 2D structured grid.

    y
    │
    │
    ○───── x

4 boundaries:

             b3
      ○─────────────○
      │             │
    b4│             │ b2
      │             │
      ○─────────────○
             b1

4 vertices:

      v4            v3
      ○─────────────○
      │             │
      │             │
      │             │
      ○─────────────○
      v1            v2

"""
meshgroupnames2D = (
    :domain_1   => :domain_1,
    :boundary_1 => :boundary_bottom,
    :boundary_2 => :boundary_right,
    :boundary_3 => :boundary_top,
    :boundary_4 => :boundary_left,
    :vertex_1   => :vertex_1,
    :vertex_2   => :vertex_2,
    :vertex_3   => :vertex_3,
    :vertex_4   => :vertex_4,
)

"""
Mesh group names on a 3D structured grid.

    z
    │ y
    │/
    ○───── x

6 boundaries:

                 b6
           ○─────/───────○
          /│            /│
         / │           / │
        /  │     b3   /  │
       ○─────────────○   │
       │b4 │         │ b2|
       │   ○─────────│───○
       │  /   b1     │  /
       │ /           │ /
       │/            │/
       ○─────────────○
              /
            b5

12 edges:

                 e7
          ○─────────────○
         /│            /│
     e8 / │e12     e6 / │
       /  │   e5     /  │e11
      ○─────────────○   │
      │   │     e3  │   │
      │   ○─────────│───○
    e9│  /       e10│  /
      │ /e4         │ / e2
      │/            │/
      ○─────────────○
             e1

8 vertices:

           v8            v7
           ○─────────────○
          /│            /│
         / │           / │
        /  │          /  │
    v5 ○─────────────○ v6│
       │   │         │   │
       │v4 ○─────────│───○ v3
       │  /          │  /
       │ /           │ /
       │/            │/
       ○─────────────○
       v1            v2
"""
meshgroupnames3D = (
    :domain_1    => :domain_1,
    :boundary_1  => :boundary_front,
    :boundary_2  => :boundary_right,
    :boundary_3  => :boundary_back,
    :boundary_4  => :boundary_left,
    :boundary_5  => :boundary_bottom,
    :boundary_6  => :boundary_top,
    :edge_1      => :edge_1,
    :edge_2      => :edge_2,
    :edge_3      => :edge_3,
    :edge_4      => :edge_4,
    :edge_5      => :edge_5,
    :edge_6      => :edge_6,
    :edge_7      => :edge_7,
    :edge_8      => :edge_8,
    :edge_9      => :edge_9,
    :edge_10     => :edge_10,
    :edge_11     => :edge_11,
    :edge_12     => :edge_12,
    :vertex_1    => :vertex_1,
    :vertex_2    => :vertex_2,
    :vertex_3    => :vertex_3,
    :vertex_4    => :vertex_4,
    :vertex_5    => :vertex_5,
    :vertex_6    => :vertex_6,
    :vertex_7    => :vertex_7,
    :vertex_8    => :vertex_8,
)

# ------------------------------------------------------------------------

Mesh(xvals, ::Type{T}=Edge2; kwargs...) where {T<:AbstractShape1D} = Mesh1D(xvals, T; kwargs...)
Mesh(xvals, yvals, ::Type{T}=Quad4; kwargs...) where {T<:AbstractShape2D} = Mesh2D(xvals, yvals, T; kwargs...)
Mesh(xvals, yvals, zvals, ::Type{T}=Hex8; kwargs...) where {T<:AbstractShape3D} = Mesh3D(xvals, yvals, zvals, T; kwargs...)

# ---------------------------------- 1D ----------------------------------

function Mesh1D(xvals, ::Type{T}=Edge2; groupnames=meshgroupnames1D, verbose=false) where {T <: AbstractShape1D}
    verbose && print("Generating 1D mesh ... ")
    coordinates = generate_coordinates(T, xvals, Float64)
    elements, domains = generate_elements(T, xvals)
    facets, boundaries = generate_facets(T, xvals)

    domain2boundaries = OrderedDict(first(keys(domains)) => collect(keys(boundaries)))

    mesh = Mesh1D(
        coordinates=coordinates,
        elements=elements,
        facets=facets,
        domains=domains,
        boundaries=boundaries,
        domain2boundaries=domain2boundaries,
        color2elements=getcolor2elements(elements),
        initcoordinates=deepcopy(coordinates),
    )

    renamegroups!(mesh, groupnames)
    checkgrouptags(mesh)

    verbose && println("Done!")

    return mesh
end

TestMesh1D(h=0.025, T=Edge3) = Mesh1D(0:h:2, T)

# ---------------------------------- 2D ----------------------------------

function Mesh2D(xvals, yvals, ::Type{T}=Quad4; groupnames=meshgroupnames2D, verbose=false) where {T <: AbstractShape2D}
    verbose && print("Generating 2D mesh ... ")
    coordinates = generate_coordinates(T, xvals, yvals, Float64)
    elements, domains = generate_elements(T, xvals, yvals)
    facets, boundaries = generate_facets(T, xvals, yvals)
    vertices, vertexgroups = generate_vertices(T, xvals, yvals)
    domain2boundaries   = OrderedDict(first(keys(domains)) => collect(keys(boundaries)))
    domain2vertexgroups = OrderedDict(first(keys(domains)) => collect(keys(vertexgroups)))

    mesh = Mesh2D(
        coordinates=coordinates,
        elements=elements,
        facets=facets,
        vertices=vertices,
        domains=domains,
        boundaries=boundaries,
        vertexgroups=vertexgroups,
        domain2boundaries=domain2boundaries,
        domain2vertexgroups=domain2vertexgroups,
        color2elements=getcolor2elements(elements),
        initcoordinates=deepcopy(coordinates),
    )

    renamegroups!(mesh, groupnames)
    checkgrouptags(mesh)

    verbose && println("Done!")

    return mesh
end

TestMesh2D(h=0.025, T=Quad8) = Mesh2D(0:h:2, 0:h:1, T)

# ---------------------------------- 3D ----------------------------------

function Mesh3D(xvals, yvals, zvals, ::Type{T}=Hex8; groupnames=meshgroupnames3D, verbose=false) where {T <: AbstractShape3D}
    verbose && print("Generating 3D mesh ... ")
    coordinates = generate_coordinates(T, xvals, yvals, zvals, Float64)
    elements, domains = generate_elements(T, xvals, yvals, zvals)
    facets, boundaries = generate_facets(T, xvals, yvals, zvals)
    edges, edgegroups = generate_edges(T, xvals, yvals, zvals)
    vertices, vertexgroups = generate_vertices(T, xvals, yvals, zvals)

    domain2boundaries   = OrderedDict(first(keys(domains)) => collect(keys(boundaries)))
    domain2edgegroups   = OrderedDict(first(keys(domains)) => collect(keys(edgegroups)))
    domain2vertexgroups = OrderedDict(first(keys(domains)) => collect(keys(vertexgroups)))

    mesh = Mesh3D(
        coordinates=coordinates,
        elements=elements,
        facets=facets,
        edges=edges,
        vertices=vertices,
        domains=domains,
        boundaries=boundaries,
        edgegroups=edgegroups,
        vertexgroups=vertexgroups,
        domain2boundaries=domain2boundaries,
        domain2edgegroups=domain2edgegroups,
        domain2vertexgroups=domain2vertexgroups,
        color2elements=getcolor2elements(elements),
        initcoordinates=deepcopy(coordinates),
    )

    renamegroups!(mesh, groupnames)
    checkgrouptags(mesh)

    verbose && println("Done!")

    return mesh
end

TestMesh3D(h=0.025, T=Hex20) = Mesh3D(0:h:2, 0:h:1, 0:h:0.5, T)

#######################
# Generate Coordinates
#######################

# ---------------------------------- 1D ----------------------------------

function generate_coordinates(::Type{Edge2}, xvals, ::Type{T}=Float64) where {T <: Number}
    coordinates = SVector{1, T}[]

    for x in xvals
        push!(coordinates, SVector{1, T}(x))
    end

    return coordinates
end

function generate_coordinates(::Type{Edge3}, xvals, ::Type{T}=Float64) where {T <: Number}
    coordinates = generate_coordinates(Edge2, xvals, T)

    for i in 1:length(xvals)-1
        x = 0.5 * (xvals[i] + xvals[i + 1])
        push!(coordinates, SVector{1, T}(x))
    end

    return coordinates
end

# ---------------------------------- 2D ----------------------------------

function generate_coordinates(::Type{<:Union{Quad4, Tri3}}, xvals, yvals, ::Type{T}=Float64) where {T <: Number}
    coordinates = SVector{2, T}[]
    for y in yvals
        for x in xvals
            push!(coordinates, SVector{2, T}(x, y))
        end
    end
    return coordinates
end

function generate_coordinates(::Type{Quad8}, xvals, yvals, ::Type{T}=Float64) where {T <: Number}
    coordinates = generate_coordinates(Quad4, xvals, yvals, T)
    for y in yvals
        for i in 1:length(xvals)-1
            x = 0.5 * (xvals[i] + xvals[i + 1])
            push!(coordinates, SVector{2, T}(x, y))
        end
    end
    for j in 1:length(yvals)-1
        for x in xvals
            y = 0.5 * (yvals[j] + yvals[j + 1])
            push!(coordinates, SVector{2, T}(x, y))
        end
    end
    return coordinates
end

function generate_coordinates(::Type{Tri6}, xvals, yvals, ::Type{T}=Float64) where {T <: Number}
    coordinates = generate_coordinates(Quad8, xvals, yvals, T)
    for j in 1:length(yvals)-1
        for i in 1:length(xvals)-1
            x = 0.5 * (xvals[i] + xvals[i + 1])
            y = 0.5 * (yvals[j] + yvals[j + 1])
            push!(coordinates, SVector{2, T}(x, y))
        end
    end
    return coordinates
end

# ---------------------------------- 3D ----------------------------------

function generate_coordinates(::Type{<:Union{Hex8, Tet4, Wedge6}}, xvals, yvals, zvals, ::Type{T}=Float64) where {T <: Number}
    coordinates = SVector{3, T}[]
    for z in zvals
        for y in yvals
            for x in xvals
                push!(coordinates, SVector{3, T}(x, y, z))
            end
        end
    end
    return coordinates
end

function generate_coordinates(::Type{Hex20}, xvals, yvals, zvals, ::Type{T}=Float64) where {T <: Number}
    coordinates = generate_coordinates(Hex8, xvals, yvals, zvals, T)
    for z in zvals
        for y in yvals
            for i in 1:length(xvals)-1
                x = 0.5 * (xvals[i] + xvals[i + 1])
                push!(coordinates, SVector{3, T}(x, y, z))
            end
        end
        for j in 1:length(yvals)-1
            for x in xvals
                y = 0.5 * (yvals[j] + yvals[j + 1])
                push!(coordinates, SVector{3, T}(x, y, z))
            end
        end
    end
    for k in 1:length(zvals)-1
        for y in yvals
            for x in xvals
                z = 0.5 * (zvals[k] + zvals[k + 1])
                push!(coordinates, SVector{3, T}(x, y, z))
            end
        end
    end
    return coordinates
end

function generate_coordinates(::Type{Tet10}, xvals, yvals, zvals, ::Type{T}=Float64) where {T <: Number}
    coordinates = generate_coordinates(Hex20, xvals, yvals, zvals, T)
    for z in zvals
        for j in 1:length(yvals)-1
            for i in 1:length(xvals)-1
                x = 0.5 * (xvals[i] + xvals[i + 1])
                y = 0.5 * (yvals[j] + yvals[j + 1])
                push!(coordinates, SVector{3, T}(x, y, z))
            end
        end
    end
    for k in 1:length(zvals)-1
        z = 0.5 * (zvals[k] + zvals[k + 1])
        for y in yvals
            for i in 1:length(xvals)-1
                x = 0.5 * (xvals[i] + xvals[i + 1])
                push!(coordinates, SVector{3, T}(x, y, z))
            end
        end
        for j in 1:length(yvals)-1
            for x in xvals
                y = 0.5 * (yvals[j] + yvals[j + 1])
                push!(coordinates, SVector{3, T}(x, y, z))
            end
        end
    end
    return coordinates
end

function generate_coordinates(::Type{Wedge15}, xvals, yvals, zvals, ::Type{T}=Float64) where {T <: Number}
    coordinates = generate_coordinates(Hex20, xvals, yvals, zvals, T)
    for z in zvals
        for j in 1:length(yvals)-1
            for i in 1:length(xvals)-1
                x = 0.5 * (xvals[i] + xvals[i + 1])
                y = 0.5 * (yvals[j] + yvals[j + 1])
                push!(coordinates, SVector{3, T}(x, y, z))
            end
        end
    end
    return coordinates
end

############################
# Grid Element Sub Division
############################

# ---------------------------------- 1D ----------------------------------
#=
    ○──────○──────○
    n1     n3     n2
=#

nsubshapes(::Type{Edge3}, ::Type{<:AbstractEdge}) = 1

function subshapeindices(::Type{Edge3}, ::Type{Edge2}, i::Integer)
    i == 1 && return SA[1, 2]
    throw(ArgumentError("No subshape Edge2 $i for shape type Edge3."))
end

function subshapeindices(::Type{Edge3}, ::Type{Edge3}, i::Integer)
    i == 1 && return SA[1, 2, 3]
    throw(ArgumentError("No subshape Edge3 $i for shape type Edge3."))
end

# ---------------------------------- 2D ----------------------------------
#=
    n4     n7     n3
    ○──────○──────○
    │             │
    ○ n8   o n9   ○ n6
    │             │
    ○──────○──────○
    n1     n5     n2
=#

nsubshapes(::Type{Quad9}, ::Type{<:AbstractQuadrilateral}) = 1
nsubshapes(::Type{Quad9}, ::Type{<:AbstractTriangle}) = 2

function subshapeindices(::Type{Quad9}, ::Type{Quad4}, i::Integer)
    i == 1 && return SA[1, 2, 3, 4]
    throw(ArgumentError("No subshape Quad4 $i for shape type Quad9."))
end

function subshapeindices(::Type{Quad9}, ::Type{Quad8}, i::Integer)
    i == 1 && return SA[1, 2, 3, 4, 5, 6, 7, 8]
    throw(ArgumentError("No subshape Quad8 $i for shape type Quad9."))
end

function subshapeindices(::Type{Quad9}, ::Type{Tri3}, i::Integer)
    i == 1 && return SA[1, 2, 3]
    i == 2 && return SA[3, 4, 1]
    throw(ArgumentError("No subshape Tri3 $i for shape type Quad9."))
end

function subshapeindices(::Type{Quad9}, ::Type{Tri6}, i::Integer)
    i == 1 && return SA[1, 2, 3, 5, 6, 9]
    i == 2 && return SA[3, 4, 1, 7, 8, 9]
    throw(ArgumentError("No subshape Tri6 $i for shape type Quad9."))
end

function rotate_nodes_z(s::Quad9)
    return Quad9(
        s.n[2], s.n[3], s.n[4], s.n[1],
        s.n[6], s.n[7], s.n[8], s.n[5],
        s.n[9],
        id=s.id,
    )
end

# ---------------------------------- 3D ----------------------------------
#=
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
=#

nsubshapes(::Type{Hex26}, ::Type{<:AbstractHexahedron}) = 1
nsubshapes(::Type{Hex26}, ::Type{<:AbstractTetrahedron}) = 5
nsubshapes(::Type{Hex26}, ::Type{<:AbstractWedge}) = 2

function subshapeindices(::Type{Hex26}, ::Type{Hex8}, i::Integer)
    i == 1 && return SA[1, 2, 3, 4, 5, 6, 7, 8]
    throw(ArgumentError("No subshape Hex8 $i for shape type Hex26."))
end

function subshapeindices(::Type{Hex26}, ::Type{Hex20}, i::Integer)
    i == 1 && return SA[1, 2, 3, 4, 5, 6, 7, 8,
        9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    throw(ArgumentError("No subshape Hex20 $i for shape type Hex26."))
end

function subshapeindices(::Type{Hex26}, ::Type{Tet4}, i::Integer)
    i == 1 && return SA[2, 3, 1, 6]
    i == 2 && return SA[7, 6, 8, 3]
    i == 3 && return SA[4, 1, 3, 8]
    i == 4 && return SA[5, 8, 6, 1]
    i == 5 && return SA[1, 8, 6, 3]
    throw(ArgumentError("No subshape Tet4 $i for shape type Hex26."))
end

function subshapeindices(::Type{Hex26}, ::Type{Tet10}, i::Integer)
    i == 1 && return SA[2, 3, 1, 6, 10, 25,  9, 18, 22, 21]
    i == 2 && return SA[7, 6, 8, 3, 14, 26, 15, 19, 22, 23]
    i == 3 && return SA[4, 1, 3, 8, 12, 25, 11, 20, 24, 23]
    i == 4 && return SA[5, 8, 6, 1, 16, 26, 13, 17, 24, 21]
    i == 5 && return SA[1, 8, 6, 3, 24, 26, 21, 25, 23, 22]
    throw(ArgumentError("No subshape Tet10 $i for shape type Hex26."))
end

function subshapeindices(::Type{Hex26}, ::Type{Wedge6}, i::Integer)
    i == 1 && return SA[1, 2, 3, 5, 6, 7]
    i == 2 && return SA[3, 4, 1, 7, 8, 5]
    throw(ArgumentError("No subshape Wedge6 $i for shape type Hex26."))
end

function subshapeindices(::Type{Hex26}, ::Type{Wedge15}, i::Integer)
    i == 1 && return SA[1, 2, 3, 5, 6, 7, 9, 10, 25, 13, 14, 26, 17, 18, 19]
    i == 2 && return SA[3, 4, 1, 7, 8, 5, 11, 12, 25, 15, 16, 26, 19, 20, 17]
    throw(ArgumentError("No subshape Wedge15 $i for shape type Hex26."))
end

function rotate_nodes_z(s::Hex26)
    return Hex26(
        s.n[2], s.n[3], s.n[4], s.n[1],
        s.n[6], s.n[7], s.n[8], s.n[5],
        s.n[10], s.n[11], s.n[12], s.n[9],
        s.n[14], s.n[15], s.n[16], s.n[13],
        s.n[18], s.n[19], s.n[20], s.n[17],
        s.n[22], s.n[23], s.n[24], s.n[21],
        s.n[25], s.n[26],
        id=s.id,
    )
end

############################
# Grid Element Constructors
############################

# ---------------------------------- 1D ----------------------------------

function GridEdge3(i::Integer, nx::Integer)
    p1 = i
    p2 = p1 + 1

    p3 = nx + i

    linear_index = i

    return Edge3(p1, p2, p3, id=linear_index)
end


# ---------------------------------- 2D ----------------------------------

function GridQuad9(i::Integer, j::Integer, nx::Integer, ny::Integer)
    # Calculations with zero-based indices
    i -= 1
    j -= 1

    nv = nx * ny
    ne = (nx - 1) * ny + nx * (ny - 1)

    p1 = j * nx + i
    p1 += 1 # Return to one-based index
    p2 = p1 + 1
    p4 = p1 + nx
    p3 = p4 + 1

    p5 = nv + j * (nx - 1) + i
    p5 += 1 # Return to one-based index
    p7 = p5 + nx - 1
    p8 = p5 + (nx - 1) * ny + j
    p6 = p8 + 1

    p9 = nv + ne + j * (nx - 1) + i
    p9 += 1 # Return to one-based index

    linear_index = (nx - 1) * j + i + 1

    return Quad9(p1, p2, p3, p4, p5, p6, p7, p8, p9, id=linear_index)
end

# ---------------------------------- 3D ----------------------------------

function GridHex26(i::Integer, j::Integer, k::Integer, nx::Integer, ny::Integer, nz::Integer)
    # Calculations with zero-based indices
    i -= 1
    j -= 1
    k -= 1

    nv = nx * ny * nz
    A = (nx - 1) * ny + nx * (ny - 1)
    ne = A * nz + nx * ny * (nz - 1)

    p1  = k * nx * ny + j * nx + i
    p1 += 1 # Return to one-based index
    p2  = p1 + 1
    p4  = p1 + nx
    p3  = p4 + 1
    p5  = p1 + nx * ny
    p6  = p5 + 1
    p8  = p5 + nx
    p7  = p8 + 1

    p9  = nv + k * A + j * (nx - 1) + i
    p9 += 1 # Return to one-based index
    p11 = p9 + nx - 1
    p12 = p9 + (nx - 1) * ny + j
    p10 = p12 + 1
    p13 = p9  + A
    p15 = p11 + A
    p16 = p12 + A
    p14 = p10 + A
    p17 = p9 + (nz - k) * A + k * nx * ny + j
    p18 = p17 + 1
    p20 = p17 + nx
    p19 = p20 + 1

    p21 = nv + ne + (nx - 1) * (ny - 1) * nz + k * A + j * (nx - 1) + i
    p21 += 1 # Return to one-based index
    p23 = p21 + nx - 1
    p24 = p21 + (nx - 1) * ny + j
    p22 = p24 + 1
    p25 = p21 + (k - nz) * (nx - 1) * (ny - 1) - k * A
    p26 = p25 + (nx - 1) * (ny - 1)

    linear_index = k * (nx - 1) * (ny - 1) + j * (nx - 1) + i + 1

    return Hex26(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16,
        p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, id=linear_index)
end

############################
# Element generators
############################

# ---------------------------------- 1D ----------------------------------

function generate_elements(::Type{T}, xvals) where {T <: AbstractShape1D}
    elements = T[]

    nx = length(xvals)

    for i in 1:nx - 1
        gridelement = GridEdge3(i, nx)
        for (m, element) in enumerate(getsubshapes(gridelement, T))
            id = i
            push!(elements, setid(element, id))
        end
    end

    domains = OrderedDict(:domain_1 => [ElementIndex(s.id) for s in elements])

    return elements, domains
end

# ---------------------------------- 2D ----------------------------------

function generate_elements(::Type{T}, xvals, yvals) where {T <: AbstractShape2D}
    elements = T[]

    nx, ny = length(xvals), length(yvals)

    for j in 1:ny - 1
        for i in 1:nx - 1
            gridelement = GridQuad9(i, j, nx, ny)
            if T<:AbstractTriangle && xor(iseven(j), iseven(i))
                gridelement = rotate_nodes_z(gridelement)
            end
            for (m, element) in enumerate(getsubshapes(gridelement, T))
                id = nsubshapes(gridelement, T) * (gridelement.id - 1) + m
                push!(elements, setid(element, id))
            end
        end
    end

    domains = OrderedDict(:domain_1 => [ElementIndex(s.id) for s in elements])

    return elements, domains
end

# ---------------------------------- 3D ----------------------------------

function generate_elements(::Type{T}, xvals, yvals, zvals) where {T <: AbstractShape3D}
    elements = T[]

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    for k in 1:nz - 1
        for j in 1:ny - 1
            for i in 1:nx - 1
                gridelement = GridHex26(i, j, k, nx, ny, nz)
                if T <: AbstractTetrahedron
                    if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
                        gridelement = rotate_nodes_z(gridelement)
                    end
                elseif T <: AbstractWedge && xor(iseven(j), iseven(i))
                    gridelement = rotate_nodes_z(gridelement)
                end
                for (m, element) in enumerate(getsubshapes(gridelement, T))
                    id = nsubshapes(gridelement, T) * (gridelement.id - 1) + m
                    push!(elements, setid(element, id))
                end
            end
        end
    end

    domains = OrderedDict(:domain_1 => [ElementIndex(s.id) for s in elements])

    return elements, domains
end

############################
# Facet generators
############################

function push_new_facet!(facets, facetindices, gridelement, T, subshape_index, localfacet_index)
    element_id = nsubshapes(gridelement, T) * (gridelement.id - 1) + subshape_index

    element = getsubshape(gridelement, T, subshape_index)
    facet = getlocalfacet(element, localfacet_index)

    facet_id = length(facets) + 1
    push!(facets, setid(facet, facet_id))
    push!(facetindices, FacetIndex(facet_id, element_id, localfacet_index))

    return facets, facetindices
end

# ---------------------------------- 1D ----------------------------------

function generate_facets(T::Type{<:AbstractEdge}, xvals)
    facets = Vertex[]
    boundaries = OrderedDict(Symbol("boundary_", i) => FacetIndex[] for i in 1:2)

    nx = length(xvals)

    push_new_facet!(facets, boundaries[:boundary_1], GridEdge3(1, nx), T, 1, 1)
    push_new_facet!(facets, boundaries[:boundary_2], GridEdge3(nx - 1, nx), T, 1, 2)

    return facets, boundaries
end

# ---------------------------------- 2D ----------------------------------

function generate_facets(T::Type{<:AbstractQuadrilateral}, xvals, yvals)
    facets = promote_type(localfacettypes(T)...)[]
    boundaries = OrderedDict(Symbol("boundary_", i) => FacetIndex[] for i in 1:4)

    nx, ny = length(xvals), length(yvals)

    j = 1
    for i in 1:nx - 1
        push_new_facet!(facets, boundaries[:boundary_1], GridQuad9(i, j, nx, ny), T, 1, 1)
    end

    i = nx - 1
    for j in 1:ny - 1
        push_new_facet!(facets, boundaries[:boundary_2], GridQuad9(i, j, nx, ny), T, 1, 2)
    end

    j = ny - 1
    for i in nx - 1:-1:1
        push_new_facet!(facets, boundaries[:boundary_3], GridQuad9(i, j, nx, ny), T, 1, 3)
    end

    i = 1
    for j in ny - 1:-1:1
        push_new_facet!(facets, boundaries[:boundary_4], GridQuad9(i, j, nx, ny), T, 1, 4)
    end

    return facets, boundaries
end

function generate_facets(T::Type{<:AbstractTriangle}, xvals, yvals)
    facets = promote_type(localfacettypes(T)...)[]
    boundaries = OrderedDict(Symbol("boundary_", i) => FacetIndex[] for i in 1:4)

    nx, ny = length(xvals), length(yvals)

    j = 1
    for i in 1:nx - 1
        gridelement = GridQuad9(i, j, nx, ny)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_facet!(facets, boundaries[:boundary_1], gridelement, T, 2, 2)
        else
            push_new_facet!(facets, boundaries[:boundary_1], gridelement, T, 1, 1)
        end
    end

    i = nx - 1
    for j in 1:ny - 1
        gridelement = GridQuad9(i, j, nx, ny)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_facet!(facets, boundaries[:boundary_2], gridelement, T, 1, 1)
        else
            push_new_facet!(facets, boundaries[:boundary_2], gridelement, T, 1, 2)
        end
    end

    j = ny - 1
    for i in nx - 1:-1:1
        gridelement = GridQuad9(i, j, nx, ny)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_facet!(facets, boundaries[:boundary_3], gridelement, T, 1, 2)
        else
            push_new_facet!(facets, boundaries[:boundary_3], gridelement, T, 2, 1)
        end
    end

    i = 1
    for j in ny - 1:-1:1
        gridelement = GridQuad9(i, j, nx, ny)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_facet!(facets, boundaries[:boundary_4], gridelement, T, 2, 1)
        else
            push_new_facet!(facets, boundaries[:boundary_4], gridelement, T, 2, 2)
        end
    end

    return facets, boundaries
end

# ---------------------------------- 3D ----------------------------------

function generate_facets(T::Type{<:AbstractHexahedron}, xvals, yvals, zvals)
    facets = promote_type(localfacettypes(T)...)[]
    boundaries = OrderedDict(Symbol("boundary_", i) => FacetIndex[] for i in 1:6)

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    j = 1
    for k in 1:nz - 1
        for i in 1:nx - 1
            push_new_facet!(facets, boundaries[:boundary_1], GridHex26(i, j, k, nx, ny, nz), T, 1, 1)
        end
    end

    i = nx - 1
    for k in 1:nz - 1
        for j in 1:ny - 1
            push_new_facet!(facets, boundaries[:boundary_2], GridHex26(i, j, k, nx, ny, nz), T, 1, 2)
        end
    end

    j = ny - 1
    for k in 1:nz - 1
        for i in 1:nx - 1
            push_new_facet!(facets, boundaries[:boundary_3], GridHex26(i, j, k, nx, ny, nz), T, 1, 3)
        end
    end

    i = 1
    for k in 1:nz - 1
        for j in 1:ny - 1
            push_new_facet!(facets, boundaries[:boundary_4], GridHex26(i, j, k, nx, ny, nz), T, 1, 4)
        end
    end

    k = 1
    for j in 1:ny - 1
        for i in 1:nx - 1
            push_new_facet!(facets, boundaries[:boundary_5], GridHex26(i, j, k, nx, ny, nz), T, 1, 5)
        end
    end

    k = nz - 1
    for j in 1:ny - 1
        for i in 1:nx - 1
            push_new_facet!(facets, boundaries[:boundary_6], GridHex26(i, j, k, nx, ny, nz), T, 1, 6)
        end
    end

    return facets, boundaries
end

function generate_facets(T::Type{<:AbstractTetrahedron}, xvals, yvals, zvals)
    facets = promote_type(localfacetypes(T)...)[]
    boundaries = OrderedDict(Symbol("boundary_", i) => FacetIndex[] for i in 1:6)

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    j = 1
    for k in 1:nz - 1
        for i in 1:nx - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_1], gridelement, T, 3, 1)
                push_new_facet!(facets, boundaries[:boundary_1], gridelement, T, 4, 1)
            else
                push_new_facet!(facets, boundaries[:boundary_1], gridelement, T, 1, 3)
                push_new_facet!(facets, boundaries[:boundary_1], gridelement, T, 4, 3)
            end
        end
    end

    i = nx - 1
    for k in 1:nz - 1
        for j in 1:ny - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_2], gridelement, T, 1, 3)
                push_new_facet!(facets, boundaries[:boundary_2], gridelement, T, 4, 3)
            else
                push_new_facet!(facets, boundaries[:boundary_2], gridelement, T, 1, 1)
                push_new_facet!(facets, boundaries[:boundary_2], gridelement, T, 2, 1)
            end
        end
    end

    j = ny - 1
    for k in 1:nz - 1
        for i in 1:nx - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_3], gridelement, T, 1, 1)
                push_new_facet!(facets, boundaries[:boundary_3], gridelement, T, 2, 1)
            else
                push_new_facet!(facets, boundaries[:boundary_3], gridelement, T, 2, 3)
                push_new_facet!(facets, boundaries[:boundary_3], gridelement, T, 3, 3)
            end
        end
    end

    i = 1
    for k in 1:nz - 1
        for j in 1:ny - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_4], gridelement, T, 2, 3)
                push_new_facet!(facets, boundaries[:boundary_4], gridelement, T, 3, 3)
            else
                push_new_facet!(facets, boundaries[:boundary_4], gridelement, T, 3, 1)
                push_new_facet!(facets, boundaries[:boundary_4], gridelement, T, 4, 1)
            end
        end
    end

    k = 1
    for j in 1:ny - 1
        for i in 1:nx - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_5], gridelement, T, 3, 4)
                push_new_facet!(facets, boundaries[:boundary_5], gridelement, T, 1, 4)
            else
                push_new_facet!(facets, boundaries[:boundary_5], gridelement, T, 1, 4)
                push_new_facet!(facets, boundaries[:boundary_5], gridelement, T, 3, 4)
            end
        end
    end

    k = nz - 1
    for j in 1:ny - 1
        for i in 1:nx - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_6], gridelement, T, 4, 4)
                push_new_facet!(facets, boundaries[:boundary_6], gridelement, T, 2, 4)
            else
                push_new_facet!(facets, boundaries[:boundary_6], gridelement, T, 2, 4)
                push_new_facet!(facets, boundaries[:boundary_6], gridelement, T, 4, 4)
            end
        end
    end

    return facets, boundaries
end

function generate_facets(T::Type{<:AbstractWedge}, xvals, yvals, zvals)
    facets = promote_type(localfacetypes(T)...)[]
    boundaries = OrderedDict(Symbol("boundary_", i) => FacetIndex[] for i in 1:6)

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    j = 1
    for k in 1:nz - 1
        for i in 1:nx - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if xor(iseven(j), iseven(i))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_1], gridelement, T, 2, 2)
            else
                push_new_facet!(facets, boundaries[:boundary_1], gridelement, T, 1, 1)
            end
        end
    end

    i = nx - 1
    for k in 1:nz - 1
        for j in 1:ny - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if xor(iseven(j), iseven(i))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_2], gridelement, T, 1, 1)
            else
                push_new_facet!(facets, boundaries[:boundary_2], gridelement, T, 1, 2)
            end
        end
    end

    j = ny - 1
    for k in 1:nz - 1
        for i in 1:nx - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if xor(iseven(j), iseven(i))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_3], gridelement, T, 1, 2)
            else
                push_new_facet!(facets, boundaries[:boundary_3], gridelement, T, 2, 1)
            end
        end
    end

    i = 1
    for k in 1:nz - 1
        for j in 1:ny - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if xor(iseven(j), iseven(i))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_4], gridelement, T, 2, 1)
            else
                push_new_facet!(facets, boundaries[:boundary_4], gridelement, T, 2, 2)
            end
        end
    end

    k = 1
    for j in 1:ny - 1
        for i in 1:nx - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if xor(iseven(j), iseven(i))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_5], gridelement, T, 2, 4)
                push_new_facet!(facets, boundaries[:boundary_5], gridelement, T, 1, 4)
            else
                push_new_facet!(facets, boundaries[:boundary_5], gridelement, T, 1, 4)
                push_new_facet!(facets, boundaries[:boundary_5], gridelement, T, 2, 4)
            end
        end
    end

    k = nz - 1
    for j in 1:ny - 1
        for i in 1:nx - 1
            gridelement = GridHex26(i, j, k, nx, ny, nz)
            if xor(iseven(j), iseven(i))
                gridelement = rotate_nodes_z(gridelement)
                push_new_facet!(facets, boundaries[:boundary_6], gridelement, T, 2, 5)
                push_new_facet!(facets, boundaries[:boundary_6], gridelement, T, 1, 5)
            else
                push_new_facet!(facets, boundaries[:boundary_6], gridelement, T, 1, 5)
                push_new_facet!(facets, boundaries[:boundary_6], gridelement, T, 2, 5)
            end
        end
    end

    return facets, boundaries
end

############################
# Edge generators
############################

function push_new_edge!(edges, edgeindices, gridelement, T, subshape_index, localedge_index)
    element_id = nsubshapes(gridelement, T) * (gridelement.id - 1) + subshape_index

    element = getsubshape(gridelement, T, subshape_index)
    edge = getlocaledge(element, localedge_index)

    edge_id = length(edges) + 1
    push!(edges, setid(edge, edge_id))
    push!(edgeindices, EdgeIndex(edge_id, element_id, localedge_index))

    return edges, edgeindices
end

# ---------------------------------- 3D ----------------------------------

function generate_edges(T::Type{<:AbstractHexahedron}, xvals, yvals, zvals)
    edges = promote_type(localedgetypes(T)...)[]
    edgegroups = OrderedDict(Symbol("edge_", i) => EdgeIndex[] for i in 1:12)

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    j, k = 1, 1
    for i in 1:nx - 1
        push_new_edge!(edges, edgegroups[:edge_1], GridHex26(i, j, k, nx, ny, nz), T, 1, 1)
    end

    i, k = nx - 1, 1
    for j in 1:ny - 1
        push_new_edge!(edges, edgegroups[:edge_2], GridHex26(i, j, k, nx, ny, nz), T, 1, 2)
    end

    j, k = ny - 1, 1
    for i in nx - 1:-1:1
        push_new_edge!(edges, edgegroups[:edge_3], GridHex26(i, j, k, nx, ny, nz), T, 1, 3)
    end

    i, k = 1, 1
    for j in ny - 1:-1:1
        push_new_edge!(edges, edgegroups[:edge_4], GridHex26(i, j, k, nx, ny, nz), T, 1, 4)
    end

    j, k = 1, nz - 1
    for i in 1:nx - 1
        push_new_edge!(edges, edgegroups[:edge_5], GridHex26(i, j, k, nx, ny, nz), T, 1, 5)
    end

    i, k = nx - 1, nz - 1
    for j in 1:ny - 1
        push_new_edge!(edges, edgegroups[:edge_6], GridHex26(i, j, k, nx, ny, nz), T, 1, 6)
    end

    j, k = ny - 1, nz - 1
    for i in nx - 1:-1:1
        push_new_edge!(edges, edgegroups[:edge_7], GridHex26(i, j, k, nx, ny, nz), T, 1, 7)
    end

    i, k = 1, nz - 1
    for j in ny - 1:-1:1
        push_new_edge!(edges, edgegroups[:edge_8], GridHex26(i, j, k, nx, ny, nz), T, 1, 8)
    end

    i, j = 1, 1
    for k in 1:nz - 1
        push_new_edge!(edges, edgegroups[:edge_9], GridHex26(i, j, k, nx, ny, nz), T, 1, 9)
    end

    i, j = nx - 1, 1
    for k in 1:nz - 1
        push_new_edge!(edges, edgegroups[:edge_10], GridHex26(i, j, k, nx, ny, nz), T, 1, 10)
    end

    i, j = nx - 1, ny - 1
    for k in 1:nz - 1
        push_new_edge!(edges, edgegroups[:edge_11], GridHex26(i, j, k, nx, ny, nz), T, 1, 11)
    end

    i, j = 1, ny - 1
    for k in 1:nz - 1
        push_new_edge!(edges, edgegroups[:edge_12], GridHex26(i, j, k, nx, ny, nz), T, 1, 12)
    end

    return edges, edgegroups
end

function generate_edges(T::Type{<:AbstractTetrahedron}, xvals, yvals, zvals)
    edges = promote_type(localedgetypes(T)...)[]
    edgegroups = OrderedDict(Symbol("edge_", i) => EdgeIndex[] for i in 1:12)

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    j, k = 1, 1
    for i in 1:nx - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_1], gridelement, T, 3, 1)
        else
            push_new_edge!(edges, edgegroups[:edge_1], gridelement, T, 1, 3)
        end
    end

    i, k = nx - 1, 1
    for j in 1:ny - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_2], gridelement, T, 1, 3)
        else
            push_new_edge!(edges, edgegroups[:edge_2], gridelement, T, 1, 1)
        end
    end

    j, k = ny - 1, 1
    for i in nx - 1:-1:1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_3], gridelement, T, 1, 1)
        else
            push_new_edge!(edges, edgegroups[:edge_3], gridelement, T, 3, 3)
        end
    end

    i, k = 1, 1
    for j in ny - 1:-1:1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_4], gridelement, T, 3, 3)
        else
            push_new_edge!(edges, edgegroups[:edge_4], gridelement, T, 3, 1)
        end
    end

    j, k = 1, nz - 1
    for i in 1:nx - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_5], gridelement, T, 4, 1)
        else
            push_new_edge!(edges, edgegroups[:edge_5], gridelement, T, 4, 3)
        end
    end

    i, k = nx - 1, nz - 1
    for j in 1:ny - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_6], gridelement, T, 4, 3)
        else
            push_new_edge!(edges, edgegroups[:edge_6], gridelement, T, 2, 1)
        end
    end

    j, k = ny - 1, nz - 1
    for i in nx - 1:-1:1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_7], gridelement, T, 2, 1)
        else
            push_new_edge!(edges, edgegroups[:edge_7], gridelement, T, 2, 3)
        end
    end

    i, k = 1, nz - 1
    for j in ny - 1:-1:1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_8], gridelement, T, 2, 3)
        else
            push_new_edge!(edges, edgegroups[:edge_8], gridelement, T, 4, 1)
        end
    end

    i, j = 1, 1
    for k in 1:nz - 1 # WARNING: Reverse edge direction!
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_9], gridelement, T, 3, 4)
        else
            push_new_edge!(edges, edgegroups[:edge_9], gridelement, T, 4, 4)
        end
    end

    i, j = nx - 1, 1
    for k in 1:nz - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_10], gridelement, T, 4, 4)
        else
            push_new_edge!(edges, edgegroups[:edge_10], gridelement, T, 1, 4)
        end
    end

    i, j = nx - 1, ny - 1
    for k in 1:nz - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_11], gridelement, T, 1, 4)
        else
            push_new_edge!(edges, edgegroups[:edge_11], gridelement, T, 2, 4)
        end
    end

    i, j = 1, ny - 1
    for k in 1:nz - 1 # WARNING: Reverse edge direction!
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_12], gridelement, T, 2, 4)
        else
            push_new_edge!(edges, edgegroups[:edge_12], gridelement, T, 3, 4)
        end
    end

    return edges, edgegroups
end

function generate_edges(T::Type{<:AbstractWedge}, xvals, yvals, zvals)
    edges = promote_type(localedgetypes(T)...)[]
    edgegroups = OrderedDict(Symbol("edge_", i) => EdgeIndex[] for i in 1:12)

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    j, k = 1, 1
    for i in 1:nx - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_1], gridelement, T, 2, 2)
        else
            push_new_edge!(edges, edgegroups[:edge_1], gridelement, T, 1, 1)
        end
    end

    i, k = nx - 1, 1
    for j in 1:ny - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_2], gridelement, T, 1, 1)
        else
            push_new_edge!(edges, edgegroups[:edge_2], gridelement, T, 1, 2)
        end
    end

    j, k = ny - 1, 1
    for i in nx - 1:-1:1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_3], gridelement, T, 1, 2)
        else
            push_new_edge!(edges, edgegroups[:edge_3], gridelement, T, 2, 1)
        end
    end

    i, k = 1, 1
    for j in ny - 1:-1:1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_4], gridelement, T, 2, 1)
        else
            push_new_edge!(edges, edgegroups[:edge_4], gridelement, T, 2, 2)
        end
    end

    j, k = 1, nz - 1
    for i in 1:nx - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_5], gridelement, T, 2, 5)
        else
            push_new_edge!(edges, edgegroups[:edge_5], gridelement, T, 1, 4)
        end
    end

    i, k = nx - 1, nz - 1
    for j in 1:ny - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_6], gridelement, T, 1, 4)
        else
            push_new_edge!(edges, edgegroups[:edge_6], gridelement, T, 1, 5)
        end
    end

    j, k = ny - 1, nz - 1
    for i in nx - 1:-1:1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_7], gridelement, T, 1, 5)
        else
            push_new_edge!(edges, edgegroups[:edge_7], gridelement, T, 2, 4)
        end
    end

    i, k = 1, nz - 1
    for j in ny - 1:-1:1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_8], gridelement, T, 2, 4)
        else
            push_new_edge!(edges, edgegroups[:edge_8], gridelement, T, 2, 5)
        end
    end

    i, j = 1, 1
    for k in 1:nz - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_9], gridelement, T, 2, 8)
        else
            push_new_edge!(edges, edgegroups[:edge_9], gridelement, T, 1, 7)
        end
    end

    i, j = nx - 1, 1
    for k in 1:nz - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_10], gridelement, T, 1, 7)
        else
            push_new_edge!(edges, edgegroups[:edge_10], gridelement, T, 1, 8)
        end
    end

    i, j = nx - 1, ny - 1
    for k in 1:nz - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_11], gridelement, T, 1, 8)
        else
            push_new_edge!(edges, edgegroups[:edge_11], gridelement, T, 2, 7)
        end
    end

    i, j = 1, ny - 1
    for k in 1:nz - 1
        gridelement = GridHex26(i, j, k, nx, ny, nz)
        if xor(iseven(j), iseven(i))
            gridelement = rotate_nodes_z(gridelement)
            push_new_edge!(edges, edgegroups[:edge_12], gridelement, T, 2, 7)
        else
            push_new_edge!(edges, edgegroups[:edge_12], gridelement, T, 2, 8)
        end
    end

    return edges, edgegroups
end

############################
# Vertex generators
############################

function push_new_vertex!(vertices, vertexindices, gridelement, T, subshape_index, localvertex_index)
    element_id = nsubshapes(gridelement, T) * (gridelement.id - 1) + subshape_index

    element = getsubshape(gridelement, T, subshape_index)
    vertex = getlocalvertex(element, localvertex_index)

    vertex_id = length(vertices) + 1
    push!(vertices, setid(vertex, vertex_id))
    push!(vertexindices, VertexIndex(vertex_id, element_id, localvertex_index))

    return vertices, vertexindices
end

# ---------------------------------- 2D ----------------------------------

function generate_vertices(T::Type{<:AbstractQuadrilateral}, xvals, yvals)
    vertices = Vertex[]
    vertexgroups = OrderedDict(Symbol("vertex_", i) => VertexIndex[] for i in 1:4)

    nx, ny, = length(xvals), length(yvals)

    i, j = 1, 1
    push_new_vertex!(vertices, vertexgroups[:vertex_1], GridQuad9(i, j, nx, ny), T, 1, 1)

    i, j = nx - 1, 1
    push_new_vertex!(vertices, vertexgroups[:vertex_2], GridQuad9(i, j, nx, ny), T, 1, 2)

    i, j = nx - 1, ny - 1
    push_new_vertex!(vertices, vertexgroups[:vertex_3], GridQuad9(i, j, nx, ny), T, 1, 3)

    i, j = 1, ny - 1
    push_new_vertex!(vertices, vertexgroups[:vertex_4], GridQuad9(i, j, nx, ny), T, 1, 4)

    return vertices, vertexgroups
end

function generate_vertices(T::Type{<:AbstractTriangle}, xvals, yvals)
    vertices = Vertex[]
    vertexgroups = OrderedDict(Symbol("vertex_", i) => VertexIndex[] for i in 1:4)

    nx, ny, = length(xvals), length(yvals)

    i, j = 1, 1
    gridelement = GridQuad9(i, j, nx, ny)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_1], gridelement, T, 2, 2)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_1], gridelement, T, 1, 1)
    end

    i, j = nx - 1, 1
    gridelement = GridQuad9(i, j, nx, ny)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_2], gridelement, T, 1, 1)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_2], gridelement, T, 1, 2)
    end

    i, j = nx - 1, ny - 1
    gridelement = GridQuad9(i, j, nx, ny)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_3], gridelement, T, 1, 2)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_3], gridelement, T, 1, 3)
    end

    i, j = 1, ny - 1
    gridelement = GridQuad9(i, j, nx, ny)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_4], gridelement, T, 1, 3)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_4], gridelement, T, 2, 2)
    end

    return vertices, vertexgroups
end

# ---------------------------------- 3D ----------------------------------

function generate_vertices(T::Type{<:AbstractHexahedron}, xvals, yvals, zvals)
    vertices = Vertex[]
    vertexgroups = OrderedDict(Symbol("vertex_", i) => VertexIndex[] for i in 1:8)

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    i, j, k = 1, 1, 1
    push_new_vertex!(vertices, vertexgroups[:vertex_1], GridHex26(i, j, k, nx, ny, nz), T, 1, 1)

    i, j, k = nx - 1, 1, 1
    push_new_vertex!(vertices, vertexgroups[:vertex_2], GridHex26(i, j, k, nx, ny, nz), T, 1, 2)

    i, j, k = nx - 1, ny - 1, 1
    push_new_vertex!(vertices, vertexgroups[:vertex_3], GridHex26(i, j, k, nx, ny, nz), T, 1, 3)

    i, j, k = 1, ny - 1, 1
    push_new_vertex!(vertices, vertexgroups[:vertex_4], GridHex26(i, j, k, nx, ny, nz), T, 1, 4)

    i, j, k = 1, 1, nz - 1
    push_new_vertex!(vertices, vertexgroups[:vertex_5], GridHex26(i, j, k, nx, ny, nz), T, 1, 5)

    i, j, k = nx - 1, 1, nz - 1
    push_new_vertex!(vertices, vertexgroups[:vertex_6], GridHex26(i, j, k, nx, ny, nz), T, 1, 6)

    i, j, k = nx - 1, ny - 1, nz - 1
    push_new_vertex!(vertices, vertexgroups[:vertex_7], GridHex26(i, j, k, nx, ny, nz), T, 1, 7)

    i, j, k = 1, ny - 1, nz - 1
    push_new_vertex!(vertices, vertexgroups[:vertex_8], GridHex26(i, j, k, nx, ny, nz), T, 1, 8)

    return vertices, vertexgroups
end


function generate_vertices(T::Type{<:AbstractTetrahedron}, xvals, yvals, zvals)
    vertices = Vertex[]
    vertexgroups = OrderedDict(Symbol("vertex_", i) => VertexIndex[] for i in 1:8)

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    i, j, k = 1, 1, 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_1], gridelement, T, 3, 1)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_1], gridelement, T, 1, 3)
    end

    i, j, k = nx - 1, 1, 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_2], gridelement, T, 1, 3)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_2], gridelement, T, 1, 1)
    end

    i, j, k = nx - 1, ny - 1, 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_3], gridelement, T, 1, 1)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_3], gridelement, T, 3, 3)
    end

    i, j, k = 1, ny - 1, 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_4], gridelement, T, 3, 3)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_4], gridelement, T, 3, 1)
    end

    i, j, k = 1, 1, nz - 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_5], gridelement, T, 2, 3)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_5], gridelement, T, 4, 1)
    end

    i, j, k = nx - 1, 1, nz - 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_6], gridelement, T, 4, 1)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_6], gridelement, T, 4, 3)
    end

    i, j, k = nx - 1, ny - 1, nz - 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_7], gridelement, T, 4, 3)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_7], gridelement, T, 2, 1)
    end

    i, j, k = 1, ny - 1, nz - 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if (isodd(k) && xor(iseven(j), iseven(i))) || (iseven(k) && !xor(iseven(j), iseven(i)))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_8], gridelement, T, 2, 1)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_8], gridelement, T, 2, 3)
    end

    return vertices, vertexgroups
end

function generate_vertices(T::Type{<:AbstractWedge}, xvals, yvals, zvals)
    vertices = Vertex[]
    vertexgroups = OrderedDict(Symbol("vertex_", i) => VertexIndex[] for i in 1:8)

    nx, ny, nz = length(xvals), length(yvals), length(zvals)

    i, j, k = 1, 1, 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_1], gridelement, T, 2, 2)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_1], gridelement, T, 1, 1)
    end

    i, j, k = nx - 1, 1, 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_2], gridelement, T, 1, 1)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_2], gridelement, T, 1, 2)
    end

    i, j, k = nx - 1, ny - 1, 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_3], gridelement, T, 1, 2)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_3], gridelement, T, 2, 1)
    end

    i, j, k = 1, ny - 1, 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_4], gridelement, T, 2, 1)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_4], gridelement, T, 2, 2)
    end

    i, j, k = 1, 1, nz - 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_5], gridelement, T, 2, 5)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_5], gridelement, T, 1, 4)
    end

    i, j, k = nx - 1, 1, nz - 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_6], gridelement, T, 1, 4)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_6], gridelement, T, 1, 5)
    end

    i, j, k = nx - 1, ny - 1, nz - 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_7], gridelement, T, 1, 5)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_7], gridelement, T, 2, 4)
    end

    i, j, k = 1, ny - 1, nz - 1
    gridelement = GridHex26(i, j, k, nx, ny, nz)
    if xor(iseven(j), iseven(i))
        gridelement = rotate_nodes_z(gridelement)
        push_new_vertex!(vertices, vertexgroups[:vertex_8], gridelement, T, 2, 4)
    else
        push_new_vertex!(vertices, vertexgroups[:vertex_8], gridelement, T, 2, 5)
    end

    return vertices, vertexgroups
end

############################
# Refinement helpers
############################

"""
    linsteprange(start, stop, startstep, stopstep)

Returns an array of `Float64` in the range from `start` to `stop`,
where the step size is varied linearly between `startstep`
and `stopstep`.

NOTE: For simplicity, step sizes `startstep` and `stopstep` are only
approximated.

Useful for mesh refinement where element sizes shall vary smoothly.

TODO: Improve algorithm while preserving symmetry for increasing
and decreasing step sizes.

# Example:
```julia-repl
julia> x1 = linsteprange(0, 1, 0.1, 0.3)
6-element Vector{Float64}:
 0.0
 0.1218467396477868
 0.2741551642075203
 0.4645406949071871
 0.7025226082817706
 1.0

julia> x2 = linsteprange(0, 1, 0.3, 0.1)
6-element Vector{Float64}:
 0.0
 0.2974773917182294
 0.5354593050928129
 0.7258448357924797
 0.8781532603522131
 1.0

julia> diff(x1) ≈ reverse(diff(x2)) # Check step size symmetry
true
```
"""
function linsteprange(start, stop, startstep, stopstep)
    slope = (stopstep - startstep) / (stop - start)

    x_up = start
    while x_up < stop
        h = startstep + slope * (x_up - start)
        x_up = x_up + h
    end

    x_down = stop
    while x_down > start
        h = stopstep + slope * (x_down - stop)
        x_down = x_down - h
    end

    if x_up - stop < start - x_down
        x = start
        x_collect = Float64[start]
        while x < stop
            h = startstep + slope * (x - start)
            x = x + h
            push!(x_collect, x)
        end

        # Final value is too large. Fix by rescaling all steps.
        x_collect .= start .+ (x_collect .- start) .* (stop - start) ./ (x_collect[end] - start)
        x_collect[1] = start
        x_collect[end] = stop

        return x_collect
    else
        x = stop
        x_collect = Float64[stop]
        while x > start
            h = stopstep + slope * (x - stop)
            x = x - h
            pushfirst!(x_collect, x)
        end

        # First value is too small. Fix by rescaling all steps.
        x_collect .= stop .+ (x_collect .- stop) .* (stop - start) ./ (stop - first(x_collect))
        x_collect[1] = start
        x_collect[end] = stop

        return x_collect
    end
end

###########################
# Annular Sector Meshes
###########################

"""
    AnnularSectorMesh(rvals, ϕvals, ::Type{T}=Quad4; kwargs...)
    AnnularSectorMesh(rvals, ϕvals, zvals, ::Type{T}=Hex8; kwargs...)

Create a mesh for the sector of an annulus (ring; difference between
concentric circles) in 2D or 3D (extrusion of an annulus into the z-direction;
difference between two coaxial cylinders).

(In 3D, this is not actually called annulus -> spherical shell)

The center of the annulus is at the origin. In 3D, the axis of the
cylinders is the z-axis.

`rvals`: Radius values r > θ

`ϕvals`: Azimuth angle values 0 < ϕ < 2π

`zvals`: Z-values (extrusion) for 3D

https://en.wikipedia.org/wiki/Annulus_(mathematics) (2D)
https://upload.wikimedia.org/wikipedia/commons/3/36/Nabla_cylindrical2.svg (3D)

# Examples
```julia-repl
julia> mesh = AnnularSectorMesh([1, 2], 0:π/4:π, Quad8);

julia> writevtk("output/annularsector2D", mesh);

julia> mesh = AnnularSectorMesh([1, 2], 0:π/4:π, [0, 1], Tet10);

julia> writevtk("output/annularsector3D", mesh);
```
"""
AnnularSectorMesh(rvals, ϕvals, ::Type{T}=Quad4; kwargs...) where {T<:AbstractShape2D} = AnnularSectorMesh2D(rvals, ϕvals, T; kwargs...)
AnnularSectorMesh(rvals, ϕvals, zvals, ::Type{T}=Hex8; kwargs...) where {T<:AbstractShape3D} = AnnularSectorMesh3D(rvals, ϕvals, zvals, T; kwargs...)

function AnnularSectorMesh2D(rvals, ϕvals, ::Type{T}=Quad4; groupnames=meshgroupnames2D, verbose=false) where {T <: AbstractShape2D}
    @assert first(rvals) < last(rvals) "rvals must be increasing"
    @assert first(rvals) > eps() "rvals must be greater than 0"
    @assert abs(last(ϕvals) - first(ϕvals)) < 2π - eps() "ϕvals must span less than 2π"
    mesh = Mesh2D(rvals, ϕvals, T; groupnames, verbose)
    T2 = eltype(mesh.coordinates)
    for (i, p) in enumerate(mesh.coordinates)
        r, ϕ = p
        mesh.coordinates[i] = T2(r * cos(ϕ), r * sin(ϕ))
    end
    mesh.initcoordinates .= mesh.coordinates
    return mesh
end

function AnnularSectorMesh3D(rvals, ϕvals, zvals, ::Type{T}=Hex8; groupnames=meshgroupnames3D, verbose=false) where {T <: AbstractShape3D}
    @assert first(rvals) < last(rvals) "rvals must be increasing"
    @assert first(rvals) > eps() "rvals must be greater than 0"
    @assert abs(last(ϕvals) - first(ϕvals)) < 2π - eps() "ϕvals must span less than 2π"
    mesh = Mesh3D(rvals, ϕvals, zvals, T; groupnames, verbose)
    T2 = eltype(mesh.coordinates)
    for (i, p) in enumerate(mesh.coordinates)
        r, ϕ, z = p
        mesh.coordinates[i] = T2(r * cos(ϕ), r * sin(ϕ), z)
    end
    mesh.initcoordinates .= mesh.coordinates
    return mesh
end

#############################
# Spherical Shell Panel Mesh
#############################

"""
    SphericalShellPanelMesh(rvals, θvals, ϕvals, ::Type{T}=Hex8; kwargs...)

Create a mesh for a panel (curved cuboid) of a spherical shell (difference between two
concentric spheres). The center of the shell is at the origin.

The panel is a curved coboid with 6 distict faces (boundaries).

`rvals`: Radius values r > θ

`θvals`: Inclination angle values 0 < θ < π

`ϕvals`: Azimuth angle values 0 < ϕ < 2π

https://upload.wikimedia.org/wikipedia/commons/3/38/Nabla_spherical2.svg

# Examples
```julia-repl
julia> mesh = SphericalShellPanelMesh([1, 2], range(π/4, 3π/4, length=5), range(0, 3π/2, length=5), Tet10);

julia> writevtk("output/sphericalshellpanel", mesh);
```
"""
function SphericalShellPanelMesh(rvals, θvals, ϕvals, ::Type{T}=Hex8; groupnames=Dict{Symbol, Symbol}(), verbose=false) where {T <: AbstractShape3D}
    @assert first(rvals) < last(rvals) "rvals must be increasing"
    @assert first(rvals) > eps() "rvals must be greater than 0"
    @assert minimum(θvals) > eps() "θvals (inclination) must be greater than 0"
    @assert maximum(θvals) < π - eps() "θvals (inclination) must be less than π"
    @assert abs(last(ϕvals) - first(ϕvals)) < 2π - eps() "ϕvals (azimuth) must span less than 2π"
    mesh = Mesh3D(rvals, θvals, ϕvals, T; groupnames, verbose)
    T2 = eltype(mesh.coordinates)
    for (i, p) in enumerate(mesh.coordinates)
        r, θ, ϕ = p
        mesh.coordinates[i] = T2(r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
    end
    mesh.initcoordinates .= mesh.coordinates
    return mesh
end
