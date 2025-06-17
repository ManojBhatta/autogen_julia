export AbstractUnvMesh, UnvMesh

##############################
# Read UNV Data Sets
##############################

"""
    readunv164!(file::IOStream)

Ignoring dataset 164 (unit data). Not implemented.

# References
https://docs.plm.automation.siemens.com/tdoc/nx/12/nx_help#uid:xid1128419:index_advanced:xid1404601:xid1404604
http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse

# Description
    UNITS
    Record 1:       FORMAT(I10,20A1,I10)
                    Field 1      -- units code
                                    = 1 - SI: Meter (newton)
                                    = 2 - BG: Foot (pound f)
                                    = 3 - MG: Meter (kilogram f)
                                    = 4 - BA: Foot (poundal)
                                    = 5 - MM: mm (milli newton)
                                    = 6 - CM: cm (centi newton)
                                    = 7 - IN: Inch (pound f)
                                    = 8 - GM: mm (kilogram f)
                                    = 9 - US: USER_DEFINED
                                    = 10- MN: mm (newton)
                    Field 2      -- units description (used for
                                    documentation only)
                    Field 3      -- temperature mode
                                    = 1 - absolute
                                    = 2 - relative
    Record 2:       FORMAT(3D25.17)
                    Unit factors for converting universal file units to SI.
                    To convert from universal file units to SI divide by
                    the appropriate factor listed below.
                    Field 1      -- length
                    Field 2      -- force
                    Field 3      -- temperature
                    Field 4      -- temperature offset
"""
function readunv164!(file::IOStream)
    while strip(readline(file)) != "-1" end
end

"""
UNV node data set
"""
struct UnvDataset2411
    nodeLabel::Vector{Int}
    coordinateSystemNo::Vector{Int}
    displacementSystemNo::Vector{Int}
    color::Vector{Int}
    xyz::Vector{SVector{3,Float64}}
end

UnvDataset2411() = UnvDataset2411(Int[], Int[], Int[], Int[], SVector{3,Float64}[])

"""
    readunv2411!(file::IOStream)

Read UNV data set 2411 (node data) and return as an instance of `UnvDataset2411`.

# References
https://docs.plm.automation.siemens.com/tdoc/nx/12/nx_help#uid:xid1128419:index_advanced:xid1404601:xid1404604
http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse

# Desciption
    NODE DATA
    Record 1:        FORMAT(4I10)
                    Field 1       -- node label
                    Field 2       -- export coordinate system number
                    Field 3       -- displacement coordinate system number
                    Field 4       -- color
    Record 2:        FORMAT(1P3D25.16)
                    Fields 1-3    -- node coordinates in the part coordinate
                                    system

    Records 1 and 2 are repeated for each node in the model.
"""
function readunv2411!(file::IOStream)
    nodeData = UnvDataset2411()

    while !eof(file)
        line1 = readline(file)

        if strip(line1) == "-1"
            break
        end

        line2 = readline(file)

        split1 = split(line1)
        split2 = split(line2)

        # Fortran scientific notation for doubles not recognized by parse()
        x = parse(Float64, replace(split2[1], "D" => "E"))
        y = parse(Float64, replace(split2[2], "D" => "E"))
        z = parse(Float64, replace(split2[3], "D" => "E"))

        push!(nodeData.nodeLabel, parse(Int, split1[1]))
        push!(nodeData.coordinateSystemNo, parse(Int, split1[2]))
        push!(nodeData.displacementSystemNo, parse(Int, split1[3]))
        push!(nodeData.color, parse(Int, split1[4]))
        push!(nodeData.xyz, SA[x, y, z])
    end

    return nodeData
end

"""
UNV element data set
"""
struct UnvDataset2412
    elementLabel::Vector{Int}
    feDescriptorId::Vector{Int}
    color::Vector{Int}
    nodeCount::Vector{Int}
    nodeLabels::Vector{Vector{Int}}
end

UnvDataset2412() = UnvDataset2412(Int[], Int[], Int[], Int[], Vector{Int}[])

"""
    unvshapetype(::Int)

Return the shape type from the unv fe descriptor id.
"""
function unvshapetype(i::Integer)
    i == 11  && return Edge2
    i == 21  && return Edge2
    i == 22  && return Edge3
    i == 24  && return Edge3
    i == 41  && return Tri3
    i == 42  && return Tri6
    i == 44  && return Quad4
    i == 45  && return Quad8
    i == 51  && return Tri3
    i == 52  && return Tri6
    i == 54  && return Quad4
    i == 55  && return Quad8
    i == 91  && return Tri3 # Actually shell element
    i == 92  && return Tri6 # Actually shell element
    i == 94  && return Quad4 # Actually shell element
    i == 95  && return Quad8 # Actually shell element
    i == 111 && return Tet4
    i == 112 && return Wedge6
    i == 113 && return Wedge15
    i == 115 && return Hex8
    i == 116 && return Hex20
    i == 118 && return Tet10
    throw(ArgumentError("unsupported unv descriptor id: $i"))
end

"""
    unv_node_perm(shape)

Return the permutation to convert a UNV shape to the node order
defined in Shapes.jl.
"""
unv_node_perm(::Type{Vertex})  = SA[1]
unv_node_perm(::Type{Edge2})   = SA[1, 2]
unv_node_perm(::Type{Edge3})   = SA[1, 3, 2]
unv_node_perm(::Type{Tri3})    = SA[1, 2, 3]
unv_node_perm(::Type{Tri6})    = SA[1, 3, 5, 2, 4, 6]
unv_node_perm(::Type{Quad4})   = SA[1, 2, 3, 4]
unv_node_perm(::Type{Quad8})   = SA[1, 3, 5, 7, 2, 4, 6, 8]
unv_node_perm(::Type{Tet4})    = SA[1, 2, 3, 4]
unv_node_perm(::Type{Tet10})   = SA[1, 3, 5, 10, 2, 4, 6, 7, 8, 9]
unv_node_perm(::Type{Hex8})    = SA[1, 2, 3, 4, 5, 6, 7, 8]
unv_node_perm(::Type{Hex20})   = SA[1, 3, 5, 7, 13, 15, 17, 19, 2, 4, 6, 8, 14, 16, 18, 20, 9, 10, 11, 12]
unv_node_perm(::Type{Wedge6})  = SA[1, 2, 3, 4, 5, 6]
unv_node_perm(::Type{Wedge15}) = SA[1, 3, 5, 10, 12, 14, 2, 4, 6, 11, 13, 15, 7, 8, 9]

"""
    readunv2412!(file::IOStream)

Read UNV data set 2412 (element data).

# References
https://docs.plm.automation.siemens.com/tdoc/nx/12/nx_help#uid:xid1128419:index_advanced:xid1404601:xid1404604
http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse

# Description
    ELEMENTS
    Record 1:        FORMAT(6I10)
                    Field 1       -- element label
                    Field 2       -- fe descriptor id
                    Field 3       -- physical property table number
                    Field 4       -- material property table number
                    Field 5       -- color
                    Field 6       -- number of nodes on element

    Record 2:  *** FOR NON-BEAM ELEMENTS ***
                    FORMAT(8I10)
                    Fields 1-n    -- node labels defining element

    Record 2:  *** FOR BEAM ELEMENTS ONLY ***
                    FORMAT(3I10)
                    Field 1       -- beam orientation node number
                    Field 2       -- beam fore-end cross section number
                    Field 3       -- beam  aft-end cross section number

    Record 3:  *** FOR BEAM ELEMENTS ONLY ***
                    FORMAT(8I10)
                    Fields 1-n    -- node labels defining element

    Records 1 and 2 are repeated for each non-beam element in the model.
    Records 1 - 3 are repeated for each beam element in the model.

    FE Descriptor Id definitions
    ____________________________
    11  Rod
    21  Linear beam
    22  Tapered beam
    23  Curved beam
    24  Parabolic beam
    31  Straight pipe
    32  Curved pipe
    41  Plane Stress Linear Triangle
    42  Plane Stress Parabolic Triangle
    43  Plane Stress Cubic Triangle
    44  Plane Stress Linear Quadrilateral
    45  Plane Stress Parabolic Quadrilateral
    46  Plane Strain Cubic Quadrilateral
    51  Plane Strain Linear Triangle
    52  Plane Strain Parabolic Triangle
    53  Plane Strain Cubic Triangle
    54  Plane Strain Linear Quadrilateral
    55  Plane Strain Parabolic Quadrilateral
    56  Plane Strain Cubic Quadrilateral
    61  Plate Linear Triangle
    62  Plate Parabolic Triangle
    63  Plate Cubic Triangle
    64  Plate Linear Quadrilateral
    65  Plate Parabolic Quadrilateral
    66  Plate Cubic Quadrilateral
    71  Membrane Linear Quadrilateral
    72  Membrane Parabolic Triangle
    73  Membrane Cubic Triangle
    74  Membrane Linear Triangle
    75  Membrane Parabolic Quadrilateral
    76  Membrane Cubic Quadrilateral
    81  Axisymetric Solid Linear Triangle
    82  Axisymetric Solid Parabolic Triangle
    84  Axisymetric Solid Linear Quadrilateral
    85  Axisymetric Solid Parabolic Quadrilateral
    91  Thin Shell Linear Triangle
    92  Thin Shell Parabolic Triangle
    93  Thin Shell Cubic Triangle
    94  Thin Shell Linear Quadrilateral
    95  Thin Shell Parabolic Quadrilateral
    96  Thin Shell Cubic Quadrilateral
    101 Thick Shell Linear Wedge
    102 Thick Shell Parabolic Wedge
    103 Thick Shell Cubic Wedge
    104 Thick Shell Linear Brick
    105 Thick Shell Parabolic Brick
    106 Thick Shell Cubic Brick
    111 Solid Linear Tetrahedron
    112 Solid Linear Wedge
    113 Solid Parabolic Wedge
    114 Solid Cubic Wedge
    115 Solid Linear Brick
    116 Solid Parabolic Brick
    117 Solid Cubic Brick
    118 Solid Parabolic Tetrahedron
    121 Rigid Bar
    122 Rigid Element
    136 Node To Node Translational Spring
    137 Node To Node Rotational Spring
    138 Node To Ground Translational Spring
    139 Node To Ground Rotational Spring
    141 Node To Node Damper
    142 Node To Gound Damper
    151 Node To Node Gap
    152 Node To Ground Gap
    161 Lumped Mass
    171 Axisymetric Linear Shell
    172 Axisymetric Parabolic Shell
    181 Constraint
    191 Plastic Cold Runner
    192 Plastic Hot Runner
    193 Plastic Water Line
    194 Plastic Fountain
    195 Plastic Baffle
    196 Plastic Rod Heater
    201 Linear node-to-node interface
    202 Linear edge-to-edge interface
    203 Parabolic edge-to-edge interface
    204 Linear face-to-face interface
    208 Parabolic face-to-face interface
    212 Linear axisymmetric interface
    213 Parabolic axisymmetric interface
    221 Linear rigid surface
    222 Parabolic rigin surface
    231 Axisymetric linear rigid surface
    232 Axisymentric parabolic rigid surface
"""
function readunv2412!(file::IOStream)
    elementData = UnvDataset2412()

    while !eof(file)
        line1 = readline(file)

        if strip(line1) == "-1"
            break
        end

        split1 = split(line1)

        # Ignoring fields 3 and 5
        push!(elementData.elementLabel, parse(Int, split1[1]))
        push!(elementData.feDescriptorId, parse(Int, split1[2]))
        push!(elementData.color, parse(Int, split1[5]))
        push!(elementData.nodeCount, parse(Int, split1[6]))

        # "Beam elements" have an additional line, which is ignored here.
        # Identify by feDescriptorId
        if elementData.feDescriptorId[end] < 40
            readline(file)
        end

        line2 = readline(file)
        split2 = split(line2)

        nodeLabels = parse.(Int, split2)

        # Continue reading nodeLabels because UNV makes a line break if there are
        # more than 8 nodes in a line.
        while length(nodeLabels) < elementData.nodeCount[end]
            line2 = readline(file)
            split2 = split(line2)
            append!(nodeLabels, parse.(Int, split2))
        end

        push!(elementData.nodeLabels, nodeLabels)
    end

    # Sort all fields by element label
    p = sortperm(elementData.elementLabel)
    elementData.elementLabel    .= elementData.elementLabel[p]
    elementData.feDescriptorId  .= elementData.feDescriptorId[p]
    elementData.color           .= elementData.color[p]
    elementData.nodeCount       .= elementData.nodeCount[p]
    elementData.nodeLabels      .= elementData.nodeLabels[p]

    return elementData
end

"""
    readunv2420!(file::IOStream)

Ignoring dataset 2420 (coordinate system data). Not implemented.

# References
https://docs.plm.automation.siemens.com/tdoc/nx/12/nx_help#uid:xid1128419:index_advanced:xid1404601:xid1404604
http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse

# Description
    COORDINATE SYSTEMS
    Record 1:        FORMAT (1I10)
                    Field 1       -- Part UID

    Record 2:        FORMAT (40A2)
                    Field 1       -- Part Name

    Record 3:        FORMAT (3I10)
                    Field 1       -- Coordinate System Label
                    Field 2       -- Coordinate System Type
                                    = 0, Cartesian
                                    = 1, Cylindrical
                                    = 2, Spherical
                    Field 3       -- Coordinate System Color

    Record 4:        FORMAT (40A2)
                    Field 1       -- Coordinate System Name

    Record 5:        FORMAT (1P3D25.16)
                    Field 1-3     -- Transformation Matrix Row 1

    Record 6:        FORMAT (1P3D25.16)
                    Field 1-3     -- Transformation Matrix Row 2

    Record 7:        FORMAT (1P3D25.16)
                    Field 1-3     -- Transformation Matrix Row 3

    Record 8:        FORMAT (1P3D25.16)
                    Field 1-3     -- Transformation Matrix Row 4

    Records 3 thru 8 are repeated for each Coordinate System in the Part.
"""
function readunv2420!(file::IOStream)
    while strip(readline(file)) != "-1" end
end

"""
UNV group data set
"""
struct UnvDataset2467
    groupNumber::Vector{Int}
    entityCount::Vector{Int}
    groupName::Vector{String}
    entityTypeCodes::Vector{Vector{Int}}
    entityTags::Vector{Vector{Int}}
end

UnvDataset2467() = UnvDataset2467(Int[], Int[], String[], Int[], Vector{Int}[])

"""
    readunv2467!(file::IOStream)

Read UNV data set 2467 (group data).

# References
https://docs.plm.automation.siemens.com/tdoc/nx/12/nx_help#uid:xid1128419:index_advanced:xid1404601:xid1404604
http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse

# Description
    PERMANENT GROUPS
    Record 1:        FORMAT(8I10)
                    Field 1       -- group number
                    Field 2       -- active constraint set no. for group
                    Field 3       -- active restraint set no. for group
                    Field 4       -- active load set no. for group
                    Field 5       -- active dof set no. for group
                    Field 6       -- active temperature set no. for group
                    Field 7       -- active contact set no. for group
                    Field 8       -- number of entities in group

    Record 2:        FORMAT(20A2)
                    Field 1       -- group name

    Record 3-N:      FORMAT(8I10)
                    Field 1       -- entity type code
                    Field 2       -- entity tag
                    Field 3       -- entity node leaf id.
                    Field 4       -- entity component/ ham id.
                    Field 5       -- entity type code
                    Field 6       -- entity tag
                    Field 7       -- entity node leaf id.
                    Field 8       -- entity component/ ham id.

    Repeat record 3 for all entities as defined by record 1, field 8.
    Records 1 thru n are repeated for each group in the model.
    Entity node leaf id. and the component/ ham id. are zero for all
    entities except "reference point", "reference point series"
    and "coordinate system".

            Permanent group entity type codes

        Entity Type Code        Entity Description

            1                coordinate system
            2                data surface thickness
            3                force on point
            4                force on edge
    ...
            7                nodes
            8                finite elements

REMAINING DESCRIPTION OMITTED
"""
function readunv2467!(file::IOStream)
    groupData = UnvDataset2467()

    while !eof(file)
        line1 = readline(file)

        if strip(line1) == "-1"
            break
        end

        split1 = split(line1)
        split2 = split(readline(file))

        # Ignoring fields 2 - 7 in record 1
        push!(groupData.groupNumber, parse(Int, split1[1]))
        entityCount = parse(Int, split1[8])
        push!(groupData.entityCount, entityCount)
        push!(groupData.groupName, split2[1])

        # Preallocate entities array
        entityTags = Array{Int}(undef, entityCount)
        entityTypeCodes = Array{Int}(undef, entityCount)

        # Fill the entities array. Each line except the last line of this group
        # contains two entities. The last line contains one entitiy if entityCount is
        # odd and two entities if entityCount is even.
        # Loop through the lines except the last one.
        readCount = 0
        if entityCount > 1
            for i = 0:div(entityCount, 2) - 1
                split3 = split(readline(file))

                entityTags[2 * i + 1] = parse(Int, split3[2])
                entityTags[2 * i + 2] = parse(Int, split3[6])
                entityTypeCodes[2 * i + 1] = parse(Int, split3[1])
                entityTypeCodes[2 * i + 2] = parse(Int, split3[5])
                readCount += 2
            end
        end

        # Now the last line if not all entities have been read
        if readCount != entityCount
            if readCount == entityCount - 1
                split3 = split(readline(file))
                entityTags[end] = parse(Int, split3[2])
                entityTypeCodes[end] = parse(Int, split3[1])
            else
                throw(ErrorException("Not all entities in unv group we were parsed!"))
            end
        end

        push!(groupData.entityTags, entityTags)
        push!(groupData.entityTypeCodes, entityTypeCodes)
    end

    return groupData
end

##############################
# UNV Mesh Type
##############################

abstract type AbstractUnvMesh end

struct UnvMesh <: AbstractUnvMesh
    nodeData::UnvDataset2411
    elementData::UnvDataset2412
    groupData::UnvDataset2467
end

Base.ndims(m::AbstractUnvMesh) = length(nonzerodims(m))

nvertices(m::AbstractUnvMesh) = length(m.nodeData.nodeLabel)
ngroups(m::AbstractUnvMesh) = length(m.groupData.groupName)

##############################
# UNV Mesh Constructors
##############################

function UnvMesh(filename::String)
    f = open(filename, "r")
    mesh = UnvMesh(f)
    close(f)
    return mesh
end

function UnvMesh(file::IOStream)
    global nodeData, elementData, groupData
    while !eof(file)
        # The unv-file contains datasets that start and end with a -1 line.
        # The second line in each dataset is a block identifier.
        if strip(readline(file)) == "-1"
            line = readline(file)
            blockid = strip(line)

            # Each of the parsing functions must read to the end of the dataset
            # and remove the ending -1 from file
            if blockid == "164"
                readunv164!(file)
            elseif blockid == "2420"
                readunv2420!(file)
            elseif blockid == "2411"
                nodeData = readunv2411!(file)
            elseif blockid == "2412"
                elementData = readunv2412!(file)
            elseif blockid in ("2467", "2477")
                groupData = readunv2467!(file)
            else
                @warn "No parsing function defined for unv dataset $blockid !"
                while strip(readline(file)) != "-1" && !eof(file) end
            end
        end
    end

    fix_elements_labels_and_node_groups!(elementData, groupData)

    return UnvMesh(nodeData, elementData, groupData)
end

#=
Set the first element label (elementData.elementLabels[1]) to 1 and reindex
groups accordingly.

If the first element label is not 1, nodes have been labeled as elements
(bug in Gmsh.jl UNV output?), but nodes are handeled differently
(entityTypeCode == 7 and not 8 in groups).
=#
function fix_elements_labels_and_node_groups!(elementData, groupData)
    offset = first(elementData.elementLabel) - 1

    if offset == 0
        return elementData, groupData
    end

    elementData.elementLabel .-= offset

    for i in eachindex(groupData.entityTypeCodes)
        # Only check the type of the first entity in group
        first(groupData.entityTypeCodes[i]) == 8 || continue

        # Only check the tag of the first entity in group
        if first(groupData.entityTags[i]) > offset
            groupData.entityTags[i] .-= offset
        else # Node group
            groupData.entityTypeCodes[i] .= 7
        end
    end

    return elementData, groupData
end

##############################
# Misc
##############################

nonzerodims(m::AbstractUnvMesh; kwargs...) = nonzerodims(m.nodeData.xyz; kwargs...)

"""
    extract_coordinates(m::AbstractUnvMesh)
    extract_coordinates(m::AbstractUnvMesh, nonzeroindices) where N

Return a vector of coordinates of UNV mesh `m`. Without argument `nonzeroindices`,
the nonzero coordinates are detected by calling `nonzerodims(m)` and the returned
coordinate vectors can therefore have length 1, 2 or 3.
Optionally, `nonzeroindices` can be passed to specifiy the length
and which coordinate indices are used, e.g. [1, 2, 3], [1, 3]  or [1].

# Examples
```julia-repl
julia> mesh3D = UnvMesh("test/meshes/mesh3D_pipe.unv");

julia> coordinates = extract_coordinates(mesh3D);

julia> coordinates[1]
3-element SVector{3, Float64} with indices SOneTo(3):
 -134.99999999999997
    0.0
  100.0

julia> mesh2D = UnvMesh("test/meshes/mesh2D_Tri3.unv");

julia> coordinates = extract_coordinates(mesh2D);

julia> coordinates[1]
2-element SVector{2, Float64} with indices SOneTo(2):
1.0
1.0

julia> coordinates = extract_coordinates(mesh2D, [1, 2, 3]);

julia> coordinates[1]
3-element SVector{3, Float64} with indices SOneTo(3):
1.0
1.0
0.0
```
"""
extract_coordinates(m::AbstractUnvMesh, args...) = extract_coordinates(m.nodeData.xyz, args...)

function extract_coordinates(xyz::Vector{<:AbstractVector}, nzdims=nonzerodims(xyz))
    @assert 1 <= length(nzdims) <= 3
    @assert minimum(nzdims) >= 1
    @assert maximum(nzdims) <= 3

    nonzero_ind = SVector{length(nzdims), Int}(nzdims...)
    coordinates = [p[nonzero_ind] for p in xyz]

    return coordinates
end

"""
    extract(unvmesh, ShapeType) -> (shapes, shapegroups)

Create and return all `shapes` of `ShapeType` in the UNV
element data in `unvmesh` (e.g. `AbstractFace`, `AbstractCell`, `Tri4`, ...)
and the `shapegroups::OrderedDict{Symbol, ElementIndex}` in which they appear.

This function loops through the groups first to match the new element indices
with unv groups. Then checks unvisited entities and attaches them.

# Note
This only returns shapes that are indexed in the unv element data.
It does not create internal edges or faces.

# Example
```julia-repl
julia> m = UnvMesh("test/meshes/mesh3D_Tet10.unv");

julia> faces, facegroups = extract(m, AbstractFace);

julia> cells, cellgroups = extract(m, AbstractCell);

julia> triangles, trianglegroups = extract(m, Tri4);
```
"""
function extract(m::AbstractUnvMesh, ShapeType::Type{<:Union{AbstractShape1D, AbstractShape2D, AbstractShape3D}})
    n_elements = length(m.elementData.elementLabel)
    groupNames = m.groupData.groupName
    entityTypeCodes = m.groupData.entityTypeCodes
    entityTags = m.groupData.entityTags

    visited = falses(n_elements)

    shapes = Vector{ShapeType}()
    shapegroups = OrderedDict{Symbol, Vector{Int}}()
    counter = 0

    # Loop through groups first
    for k in 1:ngroups(m)
        # Only check the type of the first entity in group
        first(entityTypeCodes[k]) == 8 || continue

        groupindices = Int[]

        for I in entityTags[k]
            T = unvshapetype(m.elementData.feDescriptorId[I])
            T <: ShapeType || continue

            counter += 1
            unvnodes = m.elementData.nodeLabels[I]
            shape = T(unvnodes[unv_node_perm(T)], id=counter)

            push!(shapes, shape)
            push!(groupindices, counter)
            visited[I] = true
        end

        if length(groupindices) > 0
            shapegroups[Symbol(groupNames[k])] = groupindices
        end
    end

    # Now loop through all elements that have not been used in any groups
    for (i, I) in enumerate(m.elementData.elementLabel)
        visited[i] && continue

        T = unvshapetype(m.elementData.feDescriptorId[I])
        T <: ShapeType || continue

        counter += 1
        unvnodes = m.elementData.nodeLabels[I]
        shape = T(unvnodes[unv_node_perm(T)], id=counter)

        push!(shapes, shape)
    end

    # Use comprehension to infer the element type
    shapes_ = [x for x in shapes]

    return shapes_, shapegroups
end

"""
    extract(unvmesh, Vertex) -> (vertices, shapegroups)

Create and return all `vertices` in the UNV mesh `unvmesh`
and the `vertexgroups::OrderedDict{Symbol, Int}` in which they appear.

This function loops through the groups first to match the vertex indices
with unv groups. Then checks unvisited vertices and attaches them.

# Note
This only returns vertices that are part of a UNV group. It does not
create internal vertices.

# Example
```julia-repl
julia> m = UnvMesh("test/meshes/mesh3D_Tet10.unv");

julia> vertices, vertexgroups = extract(m, AbstractVertex);
```
"""
function extract(mesh::AbstractUnvMesh, ::Type{<:AbstractShape0D})
    groupNames = mesh.groupData.groupName
    entityTypeCodes = mesh.groupData.entityTypeCodes
    entityTags = mesh.groupData.entityTags

    vertices = Vector{Vertex}()
    vertexgroups = OrderedDict{Symbol, Vector{Int}}()
    counter = 0

    # First, loop through node groups
    for k in 1:ngroups(mesh)
        # Only check the type of the first entity in group
        first(entityTypeCodes[k]) == 7 || continue

        groupindices = Int[]

        for I in entityTags[k]
            counter += 1
            vertex = Vertex(I, id=counter)

            push!(vertices, vertex)
            push!(groupindices, counter)
        end

        if length(groupindices) > 0
            vertexgroups[Symbol(groupNames[k])] = groupindices
        end
    end

    return vertices, vertexgroups
end

"""
    create_edge_elements(unvmesh) -> edges, edgegroups

Create edge elements and groups from a 1D `unvmesh`.
"""
function create_edge_elements(mesh::AbstractUnvMesh)
    elements, groups = extract(mesh, AbstractEdge)
    elementgroups = OrderedDict(name => ElementIndex.(indices) for (name, indices) in groups)
    return elements, elementgroups
end

"""
    create_face_elements(unvmesh) -> faces, facegroups

Create face elements and groups from a 2D `unvmesh`.
"""
function create_face_elements(mesh::AbstractUnvMesh)
    elements, groups = extract(mesh, AbstractFace)
    elementgroups = OrderedDict(name => ElementIndex.(indices) for (name, indices) in groups)
    return elements, elementgroups
end

"""
    create_cell_elements(unvmesh) -> cells, cellgroups

Create cell elements and groups from a 3D `unvmesh`.
"""
function create_cell_elements(mesh::AbstractUnvMesh)
    elements, groups = extract(mesh, AbstractCell)
    elementgroups = OrderedDict(name => ElementIndex.(indices) for (name, indices) in groups)
    return elements, elementgroups
end

"""
    create_group_vertices(unvmesh, elements) = vertices, element2vertices, vertexgroups

Create vertices from unv groups.

Note: Only creates those vertices contained in unv groups, not
any internal vertices.

```julia-repl
julia> unvmesh = UnvMesh("test/meshes/mesh3D_pipe.unv");

julia> elements = Felder.extract(unvmesh, AbstractCell)[1];

julia> Felder.create_group_vertices(unvmesh, elements)
```
"""
function create_group_vertices(mesh::AbstractUnvMesh, elements::Vector{<:AbstractShape})
    vertices = Vector{Vertex}()
    append_group_vertices!(vertices, mesh, elements)
end

function append_group_vertices!(vertices::Vector{<:AbstractShape0D},
                                mesh::AbstractUnvMesh, elements::Vector{<:AbstractShape1D})
    _append_group_shapes!(vertices, mesh, elements, FacetIndex)
end

function append_group_vertices!(vertices::Vector{<:AbstractShape0D},
                                mesh::AbstractUnvMesh, elements::Vector{<:Union{AbstractShape2D, AbstractShape3D}})
    _append_group_shapes!(vertices, mesh, elements, VertexIndex)
end

"""
    create_group_edges(unvmesh, elements) -> edges, element2edges, edgegroups

Create edges from unv groups.

Note: Only creates those edges contained in unv groups, not
any internal edges.

```julia-repl
julia> unvmesh = UnvMesh("test/meshes/mesh3D_pipe.unv");

julia> elements = Felder.extract(unvmesh, AbstractCell)[1];

julia> Felder.create_group_edges(unvmesh, elements)
```
"""
function create_group_edges(mesh::AbstractUnvMesh, elements::Vector{<:Union{AbstractShape2D, AbstractShape3D}})
    if eltype(elements) <: Union{Tri3, Quad4, Tet4, Hex8, Wedge6}
        T = Edge2
    elseif eltype(elements) <: Union{Tri6, Quad8, Tet10, Hex20, Wedge15}
        T = Edge3
    else
        T = AbstractEdge
    end
    edges = Vector{T}()
    append_group_edges!(edges, mesh, elements)
end

function append_group_edges!(edges::Vector{<:AbstractShape1D},
                             mesh::AbstractUnvMesh, elements::Vector{<:AbstractShape2D})
    _append_group_shapes!(edges, mesh, elements, FacetIndex)
end

function append_group_edges!(edges::Vector{<:AbstractShape1D},
                             mesh::AbstractUnvMesh, elements::Vector{<:AbstractShape3D})
    _append_group_shapes!(edges, mesh, elements, EdgeIndex)
end

"""
    create_group_faces(unvmesh, elements) -> faces, element2faces, facegroups

Create faces from unv groups.

Note: Only creates those faces contained in unv groups, not
any internal faces.

```julia-repl
julia> unvmesh = UnvMesh("test/meshes/mesh3D_pipe.unv");

julia> elements = Felder.extract(unvmesh, AbstractCell)[1];

julia> Felder.create_group_faces(unvmesh, elements)
```
"""
function create_group_faces(mesh::AbstractUnvMesh, elements::Vector{<:AbstractShape3D})
    if eltype(elements) <: Tet4
        T = Tri3
    elseif eltype(elements) <: Tet10
        T = Tri6
    elseif eltype(elements) <: Hex8
        T = Quad4
    elseif eltype(elements) <: Hex20
        T = Quad8
    else
        T = AbstractFace
    end
    faces = Vector{T}()
    append_group_faces!(faces, mesh, elements)
end

function append_group_faces!(faces::Vector{<:AbstractShape2D},
                              mesh::AbstractUnvMesh, elements::Vector{<:AbstractShape3D})
    _append_group_shapes!(faces, mesh, elements, FacetIndex)
end

function _append_group_shapes!(shapes::Vector{<:AbstractShape},
                               mesh::AbstractUnvMesh, elements::Vector{<:AbstractShape}, IndexType::Type)
    _getlocal(elements, i::FacetIndex) = getlocalfacet(elements[i.elementid], i.ilocal)
    _getlocal(elements, i::EdgeIndex) = getlocaledge(elements[i.elementid], i.ilocal)
    _getlocal(elements, i::VertexIndex) = getlocalvertex(elements[i.elementid], i.ilocal)

    node2elements = invert_map(getvertexnodes, elements)
    shapegroups = OrderedDict{Symbol, Vector{IndexType}}()
    counter = 0

    groupshapes, indexgroups = extract(mesh, eltype(shapes))

    for (groupname, indexgroup) in indexgroups
        groupindices = IndexType[]
        for i in indexgroup
            groupshape = groupshapes[i]
            i_element, i_local = find_a_parent(groupshape, elements, node2elements)

            counter += 1
            index = IndexType(counter, i_element, i_local)
            shape = _getlocal(elements, index)
            push!(shapes, setid(shape, counter))
            push!(groupindices, index)
        end
        shapegroups[groupname] = groupindices
    end

    return shapes, shapegroups
end

##############################
# UNV Mesh Constructors
##############################

function Mesh(unvmesh::AbstractUnvMesh; verbose=false, kwargs...)
    verbose && print("Reading UNV-mesh ... ")
    nzdims = nonzerodims(unvmesh)
    if length(nzdims) == 1
        m = Mesh1D(unvmesh, nzdims; kwargs...)
    elseif length(nzdims) == 2
        m = Mesh2D(unvmesh, nzdims; kwargs...)
    else
        m = Mesh3D(unvmesh, nzdims; kwargs...)
    end
    verbose && println("Done!")
    return m
end

function Mesh1D(mesh::AbstractUnvMesh, nzdims=nonzerodims(mesh))
    @assert length(nzdims) == 1

    coordinates = extract_coordinates(mesh, nzdims)
    elements, domains = create_edge_elements(mesh)

    facets, boundaries = create_group_vertices(mesh, elements)

    domain2boundaries = findsubgroups(domains, boundaries)

    Mesh1D(
        coordinates=coordinates,
        elements=elements,
        facets=facets,
        domains=domains,
        boundaries=boundaries,
        domain2boundaries=domain2boundaries,
        color2elements=getcolor2elements(elements),
        initcoordinates=deepcopy(coordinates),
    )
end

function Mesh2D(mesh::AbstractUnvMesh, nzdims=nonzerodims(mesh))
    @assert length(nzdims) == 2

    coordinates = extract_coordinates(mesh, nzdims)
    elements, domains = create_face_elements(mesh)

    vertices, vertexgroups = create_group_vertices(mesh, elements)
    facets, boundaries     = create_group_edges(mesh, elements)

    domain2vertexgroups = findsubgroups(domains, vertexgroups)
    domain2boundaries   = findsubgroups(domains, boundaries)

    Mesh2D(
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
end

function Mesh3D(mesh::AbstractUnvMesh, nzdims=nonzerodims(mesh))
    @assert length(nzdims) == 3

    coordinates = extract_coordinates(mesh, nzdims)
    elements, domains = create_cell_elements(mesh)

    vertices, vertexgroups = create_group_vertices(mesh, elements)
    edges, edgegroups      = create_group_edges(mesh, elements)
    facets, boundaries     = create_group_faces(mesh, elements)

    domain2vertexgroups = findsubgroups(domains, vertexgroups)
    domain2edgegroups   = findsubgroups(domains, edgegroups)
    domain2boundaries   = findsubgroups(domains, boundaries)

    Mesh3D(
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
end
