##############################
# Tetrahedron
##############################

center_refcoordinates(::Type{<:AbstractTetrahedron}) = SA[1/3, 1/3, 1/3]

function edge2cellcoordinates(::Type{<:AbstractTetrahedron}, i::Integer, ξ::AbstractVector)
    i == 1 && return SA[ ξ[1] + 0.5, 0.0, 0.0]
    i == 2 && return SA[-ξ[1] + 0.5, ξ[1] + 0.5, 0.0]
    i == 3 && return SA[0.0, -ξ[1] + 0.5, 0.0]
    i == 4 && return SA[0.0, 0.0, ξ[1] + 0.5]
    i == 5 && return SA[-ξ[1] + 0.5, 0.0, ξ[1] + 0.5]
    i == 6 && return SA[0.0, -ξ[1] + 0.5, ξ[1] + 0.5]
    throw(ArgumentError("No local edge $i in shape type AbstractTetrahedron"))
end

function face2cellcoordinates(::Type{<:AbstractTetrahedron}, i::Integer, ξ::AbstractVector)
    i == 1 && return SA[ξ[1], 0.0, ξ[2]]
    i == 2 && return SA[1.0 - ξ[1] - ξ[2], ξ[1], ξ[2]]
    i == 3 && return SA[0.0, 1 - ξ[1] - ξ[2], ξ[2]]
    i == 4 && return SA[ξ[2], ξ[1], 0.0]
    throw(ArgumentError("No local edge $i in shape type AbstractTetrahedron"))
end

function vertexdofs(::Type{<:AbstractTetrahedron}, ::Lagrange, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    throw(ArgumentError("No local vertex $i in AbstractTetrahedron"))
end

function refcoordinates_grid(::Type{<:AbstractTetrahedron}, n)
    # Length = sum(i -> sum(1:i), 1:n)
    return (SA[ξ, η, ζ] for (i, ζ) in enumerate(LinRange(0.0, 1.0, n))
        for (j, η) in enumerate(LinRange(0.0, 1.0 - ζ, n - i + 1))
            for ξ in LinRange(0.0, 1.0 - ζ - η, n - i - j + 2))
end

# --------------------------- Polynomial Order 1 ---------------------------
#=
     (0,0,1)
     N4 ○
        │╲⋱
        │ ╲ ⋱
        │  ╲  ⋱
    (0,1,0) ○ N3⋱
        │  /   ⋱  ⋱
        │ /      ⋱  ⋱
        │/          ⋱ ⋱
     N1 ○───────────────○ N2
     (0,0,0)         (1,0,0)
=#

ndofs(       ::Type{<:AbstractTetrahedron}, ::Lagrange{1}) = 4
ndofs_vertex(::Type{<:AbstractTetrahedron}, ::Lagrange{1}) = (1, 1, 1, 1)
ndofs_edge(  ::Type{<:AbstractTetrahedron}, ::Lagrange{1}) = (0, 0, 0, 0, 0, 0)
ndofs_face(  ::Type{<:AbstractTetrahedron}, ::Lagrange{1}) = (0, 0, 0, 0)
ndofs_cell(  ::Type{<:AbstractTetrahedron}, ::Lagrange{1}) = (0,)

function dof_refcoordinates(::Type{<:AbstractTetrahedron}, ::Lagrange{1}, i::Integer)
    # Node dofs
    i == 1 && return SA[0.0, 0.0, 0.0]
    i == 2 && return SA[1.0, 0.0, 0.0]
    i == 3 && return SA[0.0, 1.0, 0.0]
    i == 4 && return SA[0.0, 0.0, 1.0]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractTetrahedron"))
end

function shapefunc_N(::Type{<:AbstractTetrahedron}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    λ1, λ2, λ3 = λ
    λ4 = 1 - λ1 - λ2 - λ3
    # Node dofs
    i == 1 && return λ4
    i == 2 && return λ1
    i == 3 && return λ2
    i == 4 && return λ3
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractTetrahedron"))
end

function shapefunc_∇N(::Type{<:AbstractTetrahedron}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    i == 1 && return SA[-1.0,
                        -1.0,
                        -1.0]
    i == 2 && return SA[ 1.0,
                         0.0,
                         0.0]
    i == 3 && return SA[ 0.0,
                         1.0,
                         0.0]
    i == 4 && return SA[ 0.0,
                         0.0,
                         1.0]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractTetrahedron"))
end

function edgedofs(::Type{<:AbstractTetrahedron}, ::Lagrange{1}, i::Integer)
    i == 1 && return SA[1, 2]
    i == 2 && return SA[2, 3]
    i == 3 && return SA[3, 1]
    i == 4 && return SA[1, 4]
    i == 5 && return SA[2, 4]
    i == 6 && return SA[3, 4]
    throw(ArgumentError("No local edge $i in AbstractTetrahedron"))
end

function facedofs(::Type{<:AbstractTetrahedron}, ::Lagrange{1}, i::Integer)
    i == 1 && return SA[1, 2, 4]
    i == 2 && return SA[2, 3, 4]
    i == 3 && return SA[3, 1, 4]
    i == 4 && return SA[1, 3, 2]
    throw(ArgumentError("No local face $i in AbstractTetrahedron"))
end

# --------------------------- Polynomial Order 2 ---------------------------

ndofs(       ::Type{<:AbstractTetrahedron}, ::Lagrange{2}) = 10
ndofs_vertex(::Type{<:AbstractTetrahedron}, ::Lagrange{2}) = (1, 1, 1, 1)
ndofs_edge(  ::Type{<:AbstractTetrahedron}, ::Lagrange{2}) = (1, 1, 1, 1, 1, 1)
ndofs_face(  ::Type{<:AbstractTetrahedron}, ::Lagrange{2}) = (0, 0, 0, 0)
ndofs_cell(  ::Type{<:AbstractTetrahedron}, ::Lagrange{2}) = (0,)

function dof_refcoordinates(::Type{<:AbstractTetrahedron}, ::Lagrange{2}, i::Integer)
    # Node dofs
    i == 1  && return SA[0.0, 0.0, 0.0]
    i == 2  && return SA[1.0, 0.0, 0.0]
    i == 3  && return SA[0.0, 1.0, 0.0]
    i == 4  && return SA[0.0, 0.0, 1.0]
    # Edge dofs
    i == 5  && return SA[0.5, 0.0, 0.0]
    i == 6  && return SA[0.5, 0.5, 0.0]
    i == 7  && return SA[0.0, 0.5, 0.0]
    i == 8  && return SA[0.0, 0.0, 0.5]
    i == 9  && return SA[0.5, 0.0, 0.5]
    i == 10 && return SA[0.0, 0.5, 0.5]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractTetrahedron"))
end

function shapefunc_N(::Type{<:AbstractTetrahedron}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    λ1, λ2, λ3 = λ
    λ4 = 1 - λ1 - λ2 - λ3
    # Node dofs
    i == 1  && return λ4 * (2.0 * λ4 - 1)
    i == 2  && return λ1 * (2.0 * λ1 - 1)
    i == 3  && return λ2 * (2.0 * λ2 - 1)
    i == 4  && return λ3 * (2.0 * λ3 - 1)
    # Edge dofs
    i == 5  && return 4.0 * λ1 * λ4
    i == 6  && return 4.0 * λ1 * λ2
    i == 7  && return 4.0 * λ2 * λ4
    i == 8  && return 4.0 * λ3 * λ4
    i == 9  && return 4.0 * λ3 * λ1
    i == 10 && return 4.0 * λ2 * λ3
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractTetrahedron"))
end

function shapefunc_∇N(::Type{<:AbstractTetrahedron}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    λ1, λ2, λ3 = λ
    λ4 = 1 - λ1 - λ2 - λ3
    # Node dofs
    i == 1  && return SA[-4.0 * λ4 + 1,
                         -4.0 * λ4 + 1,
                         -4.0 * λ4 + 1]
    i == 2  && return SA[ 4.0 * λ1 - 1,
                          0.0,
                          0.0]
    i == 3  && return SA[ 0.0,
                          4.0 * λ2 - 1,
                          0.0]
    i == 4  && return SA[ 0.0,
                          0.0,
                          4.0 * λ3 - 1]
    # Edge dofs
    i == 5  && return SA[ 4.0 * (λ4 - λ1),
                         -4.0 * λ1,
                         -4.0 * λ1]
    i == 6  && return SA[ 4.0 * λ2,
                          4.0 * λ1,
                          0.0]
    i == 7  && return SA[-4.0 * λ2,
                          4.0 * (λ4 - λ2),
                         -4.0 * λ2]
    i == 8  && return SA[-4.0 * λ3,
                         -4.0 * λ3,
                         4.0 * (λ4 - λ3)]
    i == 9  && return SA[ 4.0 * λ3,
                          0.0,
                          4.0 * λ1]
    i == 10 && return SA[ 0.0,
                          4.0 * λ3,
                          4.0 * λ2]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractTetrahedron"))
end

function edgedofs(::Type{<:AbstractTetrahedron}, ::Lagrange{2}, i::Integer)
    i == 1 && return SA[1, 2, 5]
    i == 2 && return SA[2, 3, 6]
    i == 3 && return SA[3, 1, 7]
    i == 4 && return SA[1, 4, 8]
    i == 5 && return SA[2, 4, 9]
    i == 6 && return SA[3, 4, 10]
    throw(ArgumentError("No local edge $i in AbstractTetrahedron"))
end

function facedofs(::Type{<:AbstractTetrahedron}, ::Lagrange{2}, i::Integer)
    i == 1 && return SA[1, 2, 4, 5,  9,  8]
    i == 2 && return SA[2, 3, 4, 6, 10,  9]
    i == 3 && return SA[3, 1, 4, 7,  8, 10]
    i == 4 && return SA[1, 3, 2, 7,  6,  5]
    throw(ArgumentError("No local face $i in AbstractTetrahedron"))
end

# --------------------------- Polynomial Order 3 ---------------------------

ndofs(       ::Type{<:AbstractTetrahedron}, ::Lagrange{3}) = 20
ndofs_vertex(::Type{<:AbstractTetrahedron}, ::Lagrange{3}) = (1, 1, 1, 1)
ndofs_edge(  ::Type{<:AbstractTetrahedron}, ::Lagrange{3}) = (2, 2, 2, 2, 2, 2)
ndofs_face(  ::Type{<:AbstractTetrahedron}, ::Lagrange{3}) = (1, 1, 1, 1)
ndofs_cell(  ::Type{<:AbstractTetrahedron}, ::Lagrange{3}) = (0,)

function dof_refcoordinates(::Type{<:AbstractTetrahedron}, ::Lagrange{3}, i::Integer)
    # Node dofs
    i == 1  && return SA[0.0, 0.0, 0.0]
    i == 2  && return SA[1.0, 0.0, 0.0]
    i == 3  && return SA[0.0, 1.0, 0.0]
    i == 4  && return SA[0.0, 0.0, 1.0]
    # Edge dofs
    i == 5  && return SA[1/3, 0.0, 0.0]
    i == 6  && return SA[2/3, 0.0, 0.0]
    i == 7  && return SA[2/3, 1/3, 0.0]
    i == 8  && return SA[1/3, 2/3, 0.0]
    i == 9  && return SA[0.0, 2/3, 0.0]
    i == 10 && return SA[0.0, 1/3, 0.0]
    i == 11 && return SA[0.0, 0.0, 1/3]
    i == 12 && return SA[0.0, 0.0, 2/3]
    i == 13 && return SA[2/3, 0.0, 1/3]
    i == 14 && return SA[1/3, 0.0, 2/3]
    i == 15 && return SA[0.0, 2/3, 1/3]
    i == 16 && return SA[0.0, 1/3, 2/3]
    # Face dofs
    i == 17 && return SA[1/3, 0.0, 1/3]
    i == 18 && return SA[1/3, 1/3, 1/3]
    i == 19 && return SA[0.0, 1/3, 1/3]
    i == 20 && return SA[1/3, 1/3, 0.0]
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractTetrahedron"))
end

function shapefunc_N(::Type{<:AbstractTetrahedron}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    λ1, λ2, λ3 = λ
    λ4 = 1 - λ1 - λ2 - λ3
    # Node dofs
    i == 1  && return 0.5 * λ4 * (3 * λ4 - 1) * (3 * λ4 - 2)
    i == 2  && return 0.5 * λ1 * (3 * λ1 - 1) * (3 * λ1 - 2)
    i == 3  && return 0.5 * λ2 * (3 * λ2 - 1) * (3 * λ2 - 2)
    i == 4  && return 0.5 * λ3 * (3 * λ3 - 1) * (3 * λ3 - 2)
    # Edge dofs
    i == 5  && return 4.5 * λ1 * λ4 * (3 * λ4 - 1)
    i == 6  && return 4.5 * λ1 * λ4 * (3 * λ1 - 1)
    i == 7  && return 4.5 * λ1 * λ2 * (3 * λ1 - 1)
    i == 8  && return 4.5 * λ1 * λ2 * (3 * λ2 - 1)
    i == 9  && return 4.5 * λ2 * λ4 * (3 * λ2 - 1)
    i == 10 && return 4.5 * λ2 * λ4 * (3 * λ4 - 1)
    i == 11 && return 4.5 * λ3 * λ4 * (3 * λ4 - 1)
    i == 12 && return 4.5 * λ3 * λ4 * (3 * λ3 - 1)
    i == 13 && return 4.5 * λ3 * λ1 * (3 * λ1 - 1)
    i == 14 && return 4.5 * λ3 * λ1 * (3 * λ3 - 1)
    i == 15 && return 4.5 * λ2 * λ3 * (3 * λ2 - 1)
    i == 16 && return 4.5 * λ2 * λ3 * (3 * λ3 - 1)
    # Face dofs
    i == 17 && return 27.0 * λ3 * λ4 * λ1
    i == 18 && return 27.0 * λ1 * λ2 * λ3
    i == 19 && return 27.0 * λ2 * λ3 * λ4
    i == 20 && return 27.0 * λ4 * λ1 * λ2
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractTetrahedron"))
end

function shapefunc_∇N(::Type{<:AbstractTetrahedron}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    λ1, λ2, λ3 = λ
    λ4 = 1 - λ1 - λ2 - λ3
    # Node dofs
    i == 1  && return SA[-λ4 * (13.5 * λ4 - 9) - 1,
                         -λ4 * (13.5 * λ4 - 9) - 1,
                         -λ4 * (13.5 * λ4 - 9) - 1]
    i == 2  && return SA[λ1 * (13.5 * λ1 - 9) + 1,
                         0.0,
                         0.0]
    i == 3  && return SA[0.0,
                         λ2 * (13.5 * λ2 - 9) + 1,
                         0.0]
    i == 4  && return SA[0.0,
                         0.0,
                         λ3 * (13.5 * λ3 - 9) + 1]
    # Edge dofs
    i == 5  && return SA[  4.5 * (λ4 * ( 3 * λ4 - 6 * λ1 - 1) + λ1),
                         -27.0 * λ1 * λ4 + 4.5 * λ1,
                         -27.0 * λ1 * λ4 + 4.5 * λ1]
    i == 6  && return SA[  4.5 * (λ1 * (-3 * λ1 + 6 * λ4 + 1) - λ4),
                          -4.5 * λ1 * (3 * λ1 - 1),
                          -4.5 * λ1 * (3 * λ1 - 1)]
    i == 7  && return SA[ 27.0 * λ1 * λ2 - 4.5 * λ2,
                           4.5 * λ1 * (3 * λ1 - 1),
                           0.0]
    i == 8  && return SA[  4.5 * λ2 * (3 * λ2 - 1),
                          27.0 * λ1 * λ2 - 4.5 * λ1,
                           0.0]
    i == 9  && return SA[ -4.5 * λ2 * (3 * λ2 - 1),
                           4.5 * (λ2 * (-3 * λ2 + 6 * λ4 + 1) - λ4),
                          -4.5 * λ2 * (3 * λ2 - 1)]
    i == 10 && return SA[-27.0 * λ2 * λ4 + 4.5 * λ2,
                           4.5 * (λ4 * ( 3 * λ4 - 6 * λ2 - 1) + λ2),
                         -27.0 * λ2 * λ4 + 4.5 * λ2]
    i == 11 && return SA[-27.0 * λ3 * λ4 + 4.5 * λ3,
                         -27.0 * λ3 * λ4 + 4.5 * λ3,
                           4.5 * (λ4 * ( 3 * λ4 - 6 * λ3 - 1) + λ3)]
    i == 12 && return SA[ -4.5 * λ3 * (3 * λ3 - 1),
                          -4.5 * λ3 * (3 * λ3 - 1),
                           4.5 * (λ3 * (-3 * λ3 + 6 * λ4 + 1) - λ4)]
    i == 13 && return SA[ 27.0 * λ3 * λ1 - 4.5 * λ3,
                           0.0,
                           4.5 * λ1 * (3 * λ1 - 1)]
    i == 14 && return SA[  4.5 * λ3 * (3 * λ3 - 1),
                           0.0,
                          27.0 * λ3 * λ1 - 4.5 * λ1]
    i == 15 && return SA[  0.0,
                          27.0 * λ2 * λ3 - 4.5 * λ3,
                           4.5 * λ2 * (3 * λ2 - 1)]
    i == 16 && return SA[  0.0,
                           4.5 * λ3 * (3 * λ3 - 1),
                          27.0 * λ2 * λ3 - 4.5 * λ2]
    # Face dofs
    i == 17 && return SA[ 27.0 * λ3 * (λ4 - λ1),
                         -27.0 * λ3 * λ1,
                          27.0 * λ1 * (λ4 - λ3)]
    i == 18 && return SA[ 27.0 * λ2 * λ3,
                          27.0 * λ1 * λ3,
                          27.0 * λ1 * λ2]
    i == 19 && return SA[-27.0 * λ2 * λ3,
                          27.0 * λ3 * (λ4 - λ2),
                          27.0 * λ2 * (λ4 - λ3)]
    i == 20 && return SA[ 27.0 * λ2 * (λ4 - λ1),
                          27.0 * λ1 * (λ4 - λ2),
                         -27.0 * λ1 * λ2]

    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractTetrahedron"))
end

function edgedofs(::Type{<:AbstractTetrahedron}, ::Lagrange{3}, i::Integer)
    i == 1 && return SA[1, 2, 5, 6]
    i == 2 && return SA[2, 3, 7, 8]
    i == 3 && return SA[3, 1, 9, 10]
    i == 4 && return SA[1, 4, 11, 12]
    i == 5 && return SA[2, 4, 13, 14]
    i == 6 && return SA[3, 4, 15, 16]
    throw(ArgumentError("No local edge $i in AbstractTetrahedron"))
end

function facedofs(::Type{<:AbstractTetrahedron}, ::Lagrange{3}, i::Integer)
    i == 1 && return SA[1, 2, 4, 5, 6, 13, 14, 12, 11, 17]
    i == 2 && return SA[2, 3, 4, 7, 8, 15, 16, 14, 13, 18]
    i == 3 && return SA[3, 1, 4, 9, 10, 11, 12, 16, 15, 19]
    i == 4 && return SA[1, 3, 2, 10, 9, 8, 7, 6, 5, 20]
    throw(ArgumentError("No local face $i in AbstractTetrahedron"))
end

##############################
# Hexahedron
##############################

center_refcoordinates(::Type{<:AbstractHexahedron}) = SA[0.0, 0.0, 0.0]

function edge2cellcoordinates(::Type{<:AbstractHexahedron}, i::Integer, ξ::AbstractVector)
    i ==  1 && return SA[ ξ[1], -0.5, -0.5]
    i ==  2 && return SA[ 0.5,  ξ[1], -0.5]
    i ==  3 && return SA[-ξ[1],  0.5, -0.5]
    i ==  4 && return SA[-0.5, -ξ[1], -0.5]
    i ==  5 && return SA[ ξ[1], -0.5,  0.5]
    i ==  6 && return SA[ 0.5,  ξ[1],  0.5]
    i ==  7 && return SA[-ξ[1],  0.5,  0.5]
    i ==  8 && return SA[-0.5, -ξ[1],  0.5]
    i ==  9 && return SA[-0.5, -0.5,  ξ[1]]
    i == 10 && return SA[ 0.5, -0.5,  ξ[1]]
    i == 11 && return SA[ 0.5,  0.5,  ξ[1]]
    i == 12 && return SA[-0.5,  0.5,  ξ[1]]
    throw(ArgumentError("No local edge $i in shape type AbstractHexahedron"))
end

function face2cellcoordinates(::Type{<:AbstractHexahedron}, i::Integer, ξ::AbstractVector)
    i == 1 && return SA[ ξ[1],  -0.5,  ξ[2]]
    i == 2 && return SA[  0.5,  ξ[1],  ξ[2]]
    i == 3 && return SA[-ξ[1],   0.5,  ξ[2]]
    i == 4 && return SA[ -0.5, -ξ[1],  ξ[2]]
    i == 5 && return SA[ ξ[2],  ξ[1], -0.5]
    i == 6 && return SA[ ξ[1],  ξ[2],  0.5]
    throw(ArgumentError("No local edge $i in shape type AbstractHexahedron"))
end

function vertexdofs(::Type{<:AbstractHexahedron}, ::Lagrange, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    i == 5 && return 5
    i == 6 && return 6
    i == 7 && return 7
    i == 8 && return 8
    throw(ArgumentError("No local vertex $i in AbstractHexahedron"))
end

function refcoordinates_grid(::Type{<:AbstractHexahedron}, n)
    # Length = n * n * n
    return (SA[ξ, η, ζ] for ζ in LinRange(-0.5, 0.5, n)
        for η in LinRange(-0.5, 0.5, n)
            for ξ in LinRange(-0.5, 0.5, n))
end

# --------------------------- Polynomial Order 1 ---------------------------
#=
Linear Hexahedron defined by 8 points

             (-½,½,½)        (½,½,½)
                8 ○─────────────○ 7
                 /│            /│
                / │           / │
             5 /  │          /  │
    (-½,-½,½) ○─────────────○ 6 │
              │   │         │   │
              │ 4 ○─────────│───○ (½,½,-½)
              │  /          │  / 3
              │ /           │ /
              │/            │/
            1 ○─────────────○ 2
         (-½,-½,-½)     (½,-½,-½)
=#

ndofs(       ::Type{<:AbstractHexahedron}, ::Lagrange{1}) = 8
ndofs_vertex(::Type{<:AbstractHexahedron}, ::Lagrange{1}) = (1, 1, 1, 1, 1, 1, 1, 1)
ndofs_edge(  ::Type{<:AbstractHexahedron}, ::Lagrange{1}) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
ndofs_face(  ::Type{<:AbstractHexahedron}, ::Lagrange{1}) = (0, 0, 0, 0, 0, 0)
ndofs_cell(  ::Type{<:AbstractHexahedron}, ::Lagrange{1}) = (0,)

function dof_refcoordinates(::Type{<:AbstractHexahedron}, ::Lagrange{1}, i::Integer)
    # Node dofs
    i == 1 && return SA[-0.5, -0.5, -0.5]
    i == 2 && return SA[ 0.5, -0.5, -0.5]
    i == 3 && return SA[ 0.5,  0.5, -0.5]
    i == 4 && return SA[-0.5,  0.5, -0.5]
    i == 5 && return SA[-0.5, -0.5,  0.5]
    i == 6 && return SA[ 0.5, -0.5,  0.5]
    i == 7 && return SA[ 0.5,  0.5,  0.5]
    i == 8 && return SA[-0.5,  0.5,  0.5]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractHexahedron"))
end

function shapefunc_N(::Type{<:AbstractHexahedron}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    x, y, z = λ
    # Node dofs
    i == 1 && return (0.5 - x) * (0.5 - y) * (0.5 - z)
    i == 2 && return (0.5 + x) * (0.5 - y) * (0.5 - z)
    i == 3 && return (0.5 + x) * (0.5 + y) * (0.5 - z)
    i == 4 && return (0.5 - x) * (0.5 + y) * (0.5 - z)
    i == 5 && return (0.5 - x) * (0.5 - y) * (0.5 + z)
    i == 6 && return (0.5 + x) * (0.5 - y) * (0.5 + z)
    i == 7 && return (0.5 + x) * (0.5 + y) * (0.5 + z)
    i == 8 && return (0.5 - x) * (0.5 + y) * (0.5 + z)
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractHexahedron"))
end

function shapefunc_∇N(::Type{<:AbstractHexahedron}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    x, y, z = λ
    # Node dofs
    i == 1  && return SA[-(0.5 - y) * (0.5 - z),
                         -(0.5 - x) * (0.5 - z),
                         -(0.5 - x) * (0.5 - y)]
    i == 2  && return SA[ (0.5 - y) * (0.5 - z),
                         -(0.5 + x) * (0.5 - z),
                         -(0.5 + x) * (0.5 - y)]
    i == 3  && return SA[ (0.5 + y) * (0.5 - z),
                          (0.5 + x) * (0.5 - z),
                         -(0.5 + x) * (0.5 + y)]
    i == 4  && return SA[-(0.5 + y) * (0.5 - z),
                          (0.5 - x) * (0.5 - z),
                         -(0.5 - x) * (0.5 + y)]
    i == 5  && return SA[-(0.5 - y) * (0.5 + z),
                         -(0.5 - x) * (0.5 + z),
                          (0.5 - x) * (0.5 - y)]
    i == 6  && return SA[ (0.5 - y) * (0.5 + z),
                         -(0.5 + x) * (0.5 + z),
                          (0.5 + x) * (0.5 - y)]
    i == 7  && return SA[ (0.5 + y) * (0.5 + z),
                          (0.5 + x) * (0.5 + z),
                          (0.5 + x) * (0.5 + y)]
    i == 8  && return SA[-(0.5 + y) * (0.5 + z),
                          (0.5 - x) * (0.5 + z),
                          (0.5 - x) * (0.5 + y)]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractHexahedron"))
end

function edgedofs(::Type{<:AbstractHexahedron}, ::Lagrange{1}, i::Integer)
    i == 1  && return SA[1, 2]
    i == 2  && return SA[2, 3]
    i == 3  && return SA[3, 4]
    i == 4  && return SA[4, 1]
    i == 5  && return SA[5, 6]
    i == 6  && return SA[6, 7]
    i == 7  && return SA[7, 8]
    i == 8  && return SA[8, 5]
    i == 9  && return SA[1, 5]
    i == 10 && return SA[2, 6]
    i == 11 && return SA[3, 7]
    i == 12 && return SA[4, 8]
    throw(ArgumentError("No local edge $i in AbstractHexahedron"))
end

function facedofs(::Type{<:AbstractHexahedron}, ::Lagrange{1}, i::Integer)
    i == 1 && return SA[1, 2, 6, 5]
    i == 2 && return SA[2, 3, 7, 6]
    i == 3 && return SA[3, 4, 8, 7]
    i == 4 && return SA[4, 1, 5, 8]
    i == 5 && return SA[1, 4, 3, 2]
    i == 6 && return SA[5, 6, 7, 8]
    throw(ArgumentError("No local face $i in AbstractHexahedron"))
end

# --------------------------- Polynomial Order 2 ---------------------------

ndofs(       ::Type{<:AbstractHexahedron}, ::Lagrange{2}) = 20
ndofs_vertex(::Type{<:AbstractHexahedron}, ::Lagrange{2}) = (1, 1, 1, 1, 1, 1, 1, 1,)
ndofs_edge(  ::Type{<:AbstractHexahedron}, ::Lagrange{2}) = (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
ndofs_face(  ::Type{<:AbstractHexahedron}, ::Lagrange{2}) = (0, 0, 0, 0, 0, 0)
ndofs_cell(  ::Type{<:AbstractHexahedron}, ::Lagrange{2}) = (0,)

function dof_refcoordinates(::Type{<:AbstractHexahedron}, ::Lagrange{2}, i::Integer)
    # Node dofs
    i == 1  && return SA[-0.5, -0.5, -0.5]
    i == 2  && return SA[ 0.5, -0.5, -0.5]
    i == 3  && return SA[ 0.5,  0.5, -0.5]
    i == 4  && return SA[-0.5,  0.5, -0.5]
    i == 5  && return SA[-0.5, -0.5,  0.5]
    i == 6  && return SA[ 0.5, -0.5,  0.5]
    i == 7  && return SA[ 0.5,  0.5,  0.5]
    i == 8  && return SA[-0.5,  0.5,  0.5]
    # Edge dofs
    i == 9  && return SA[ 0.0, -0.5, -0.5]
    i == 10 && return SA[ 0.5,  0.0, -0.5]
    i == 11 && return SA[ 0.0,  0.5, -0.5]
    i == 12 && return SA[-0.5,  0.0, -0.5]
    i == 13 && return SA[ 0.0, -0.5,  0.5]
    i == 14 && return SA[ 0.5,  0.0,  0.5]
    i == 15 && return SA[ 0.0,  0.5,  0.5]
    i == 16 && return SA[-0.5,  0.0,  0.5]
    i == 17 && return SA[-0.5, -0.5,  0.0]
    i == 18 && return SA[ 0.5, -0.5,  0.0]
    i == 19 && return SA[ 0.5,  0.5,  0.0]
    i == 20 && return SA[-0.5,  0.5,  0.0]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractHexahedron"))
end

function shapefunc_N(::Type{<:AbstractHexahedron}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    x, y, z = λ .* 2
    # Node dofs
    i == 1  && return 0.125 * (1 - x) * (1 - y) * (1 - z) * (-x - y - z - 2)
    i == 2  && return 0.125 * (1 + x) * (1 - y) * (1 - z) * ( x - y - z - 2)
    i == 3  && return 0.125 * (1 + x) * (1 + y) * (1 - z) * ( x + y - z - 2)
    i == 4  && return 0.125 * (1 - x) * (1 + y) * (1 - z) * (-x + y - z - 2)
    i == 5  && return 0.125 * (1 - x) * (1 - y) * (1 + z) * (-x - y + z - 2)
    i == 6  && return 0.125 * (1 + x) * (1 - y) * (1 + z) * ( x - y + z - 2)
    i == 7  && return 0.125 * (1 + x) * (1 + y) * (1 + z) * ( x + y + z - 2)
    i == 8  && return 0.125 * (1 - x) * (1 + y) * (1 + z) * (-x + y + z - 2)
    # Edge dofs
    i == 9  && return 0.25 * (1 - x^2) * (1 - y)   * (1 - z)
    i == 10 && return 0.25 * (1 + x)   * (1 - y^2) * (1 - z)
    i == 11 && return 0.25 * (1 - x^2) * (1 + y)   * (1 - z)
    i == 12 && return 0.25 * (1 - x)   * (1 - y^2) * (1 - z)
    i == 13 && return 0.25 * (1 - x^2) * (1 - y)   * (1 + z)
    i == 14 && return 0.25 * (1 + x)   * (1 - y^2) * (1 + z)
    i == 15 && return 0.25 * (1 - x^2) * (1 + y)   * (1 + z)
    i == 16 && return 0.25 * (1 - x)   * (1 - y^2) * (1 + z)
    i == 17 && return 0.25 * (1 - x)   * (1 - y)   * (1 - z^2)
    i == 18 && return 0.25 * (1 + x)   * (1 - y)   * (1 - z^2)
    i == 19 && return 0.25 * (1 + x)   * (1 + y)   * (1 - z^2)
    i == 20 && return 0.25 * (1 - x)   * (1 + y)   * (1 - z^2)
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractHexahedron"))
end


function shapefunc_∇N(::Type{<:AbstractHexahedron}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    x, y, z = λ .* 2
    # Node dofs
    i == 1  && return SA[ 0.25 * (1 - y) * (1 - z) * ( 2*x +   y +   z + 1),
                          0.25 * (1 - x) * (1 - z) * (   x + 2*y +   z + 1),
                          0.25 * (1 - x) * (1 - y) * (   x +   y + 2*z + 1)]
    i == 2  && return SA[ 0.25 * (1 - y) * (1 - z) * ( 2*x -   y -   z - 1),
                         -0.25 * (1 + x) * (1 - z) * (   x - 2*y -   z - 1),
                         -0.25 * (1 + x) * (1 - y) * (   x -   y - 2*z - 1)]
    i == 3  && return SA[ 0.25 * (1 + y) * (1 - z) * ( 2*x +   y -   z - 1),
                          0.25 * (1 + x) * (1 - z) * (   x + 2*y -   z - 1),
                         -0.25 * (1 + x) * (1 + y) * (   x +   y - 2*z - 1)]
    i == 4  && return SA[-0.25 * (1 + y) * (1 - z) * (-2*x +   y -   z - 1),
                          0.25 * (1 - x) * (1 - z) * (  -x + 2*y -   z - 1),
                         -0.25 * (1 - x) * (1 + y) * (  -x +   y - 2*z - 1)]
    i == 5  && return SA[-0.25 * (1 - y) * (1 + z) * (-2*x -   y +   z - 1),
                         -0.25 * (1 - x) * (1 + z) * (  -x - 2*y +   z - 1),
                          0.25 * (1 - x) * (1 - y) * (  -x -   y + 2*z - 1)]
    i == 6  && return SA[ 0.25 * (1 - y) * (1 + z) * ( 2*x -   y +   z - 1),
                         -0.25 * (1 + x) * (1 + z) * (   x - 2*y +   z - 1),
                          0.25 * (1 + x) * (1 - y) * (   x -   y + 2*z - 1)]
    i == 7  && return SA[ 0.25 * (1 + y) * (1 + z) * ( 2*x +   y +   z - 1),
                          0.25 * (1 + x) * (1 + z) * (   x + 2*y +   z - 1),
                          0.25 * (1 + x) * (1 + y) * (   x +   y + 2*z - 1)]
    i == 8  && return SA[-0.25 * (1 + y) * (1 + z) * (-2*x +   y +   z - 1),
                          0.25 * (1 - x) * (1 + z) * (  -x + 2*y +   z - 1),
                          0.25 * (1 - x) * (1 + y) * (  -x +   y + 2*z - 1)]
    # Edge dofs
    i == 9  && return SA[ -x * (1 - y) * (1 - z),
                         -0.5 * (1 - x^2) * (1 - z),
                         -0.5 * (1 - x^2) * (1 - y)]
    i == 10 && return SA[ 0.5 * (1 - y^2) * (1 - z),
                          -y * (1 + x) * (1 - z),
                         -0.5 * (1 + x) * (1 - y^2)]
    i == 11 && return SA[ -x * (1 + y) * (1 - z),
                          0.5 * (1 - x^2) * (1 - z),
                         -0.5 * (1 - x^2) * (1 + y)]
    i == 12 && return SA[-0.5 * (1 - y^2) * (1 - z),
                          -y * (1 - x) * (1 - z),
                         -0.5 * (1 - x) * (1 - y^2)]
    i == 13 && return SA[ -x * (1 - y) * (1 + z),
                         -0.5 * (1 - x^2) * (1 + z),
                          0.5 * (1 - x^2) * (1 - y)]
    i == 14 && return SA[ 0.5 * (1 - y^2) * (1 + z),
                          -y * (1 + x) * (1 + z),
                          0.5 * (1 + x) * (1 - y^2)]
    i == 15 && return SA[ -x * (1 + y) * (1 + z),
                          0.5 * (1 - x^2) * (1 + z),
                          0.5 * (1 - x^2) * (1 + y)]
    i == 16 && return SA[-0.5 * (1 - y^2) * (1 + z),
                          -y * (1 - x) * (1 + z),
                          0.5 * (1 - x) * (1 - y^2)]
    i == 17 && return SA[-0.5 * (1 - y) * (1 - z^2),
                         -0.5 * (1 - x) * (1 - z^2),
                          -z * (1 - x) * (1 - y)]
    i == 18 && return SA[ 0.5 * (1 - y) * (1 - z^2),
                         -0.5 * (1 + x) * (1 - z^2),
                          -z * (1 + x) * (1 - y)]
    i == 19 && return SA[ 0.5 * (1 + y) * (1 - z^2),
                          0.5 * (1 + x) * (1 - z^2),
                          -z * (1 + x) * (1 + y)]
    i == 20 && return SA[-0.5 * (1 + y) * (1 - z^2),
                          0.5 * (1 - x) * (1 - z^2),
                          -z * (1 - x) * (1 + y)]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractHexahedron"))
end

function edgedofs(::Type{<:AbstractHexahedron}, ::Lagrange{2}, i::Integer)
    i == 1  && return SA[1, 2,  9]
    i == 2  && return SA[2, 3, 10]
    i == 3  && return SA[3, 4, 11]
    i == 4  && return SA[4, 1, 12]
    i == 5  && return SA[5, 6, 13]
    i == 6  && return SA[6, 7, 14]
    i == 7  && return SA[7, 8, 15]
    i == 8  && return SA[8, 5, 16]
    i == 9  && return SA[1, 5, 17]
    i == 10 && return SA[2, 6, 18]
    i == 11 && return SA[3, 7, 19]
    i == 12 && return SA[4, 8, 20]
    throw(ArgumentError("No local edge $i in AbstractHexahedron"))
end

function facedofs(::Type{<:AbstractHexahedron}, ::Lagrange{2}, i::Integer)
    i == 1 && return SA[1, 2, 6, 5,  9, 18, 13, 17]
    i == 2 && return SA[2, 3, 7, 6, 10, 19, 14, 18]
    i == 3 && return SA[3, 4, 8, 7, 11, 20, 15, 19]
    i == 4 && return SA[4, 1, 5, 8, 12, 17, 16, 20]
    i == 5 && return SA[1, 4, 3, 2, 12, 11, 10,  9]
    i == 6 && return SA[5, 6, 7, 8, 13, 14, 15, 16]
    throw(ArgumentError("No local face $i in AbstractHexahedron"))
end
# --------------------------- Polynomial Order 3 ---------------------------

ndofs(       ::Type{<:AbstractHexahedron}, ::Lagrange{3}) = 32
ndofs_vertex(::Type{<:AbstractHexahedron}, ::Lagrange{3}) = (1, 1, 1, 1, 1, 1, 1, 1,)
ndofs_edge(  ::Type{<:AbstractHexahedron}, ::Lagrange{3}) = (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
ndofs_face(  ::Type{<:AbstractHexahedron}, ::Lagrange{3}) = (0, 0, 0, 0, 0, 0)
ndofs_cell(  ::Type{<:AbstractHexahedron}, ::Lagrange{3}) = (0,)

function dof_refcoordinates(::Type{<:AbstractHexahedron}, ::Lagrange{3}, i::Integer)
    # Node dofs
    i == 1  && return SA[-0.5, -0.5, -0.5]
    i == 2  && return SA[ 0.5, -0.5, -0.5]
    i == 3  && return SA[ 0.5,  0.5, -0.5]
    i == 4  && return SA[-0.5,  0.5, -0.5]
    i == 5  && return SA[-0.5, -0.5,  0.5]
    i == 6  && return SA[ 0.5, -0.5,  0.5]
    i == 7  && return SA[ 0.5,  0.5,  0.5]
    i == 8  && return SA[-0.5,  0.5,  0.5]
    # Edge dofs
    i == 9  && return SA[-1/6, -0.5, -0.5]
    i == 10 && return SA[ 1/6, -0.5, -0.5]
    i == 11 && return SA[ 0.5, -1/6, -0.5]
    i == 12 && return SA[ 0.5,  1/6, -0.5]
    i == 13 && return SA[ 1/6,  0.5, -0.5]
    i == 14 && return SA[-1/6,  0.5, -0.5]
    i == 15 && return SA[-0.5,  1/6, -0.5]
    i == 16 && return SA[-0.5, -1/6, -0.5]
    i == 17 && return SA[-1/6, -0.5,  0.5]
    i == 18 && return SA[ 1/6, -0.5,  0.5]
    i == 19 && return SA[ 0.5, -1/6,  0.5]
    i == 20 && return SA[ 0.5,  1/6,  0.5]
    i == 21 && return SA[ 1/6,  0.5,  0.5]
    i == 22 && return SA[-1/6,  0.5,  0.5]
    i == 23 && return SA[-0.5,  1/6,  0.5]
    i == 24 && return SA[-0.5, -1/6,  0.5]
    i == 25 && return SA[-0.5, -0.5, -1/6]
    i == 26 && return SA[-0.5, -0.5,  1/6]
    i == 27 && return SA[ 0.5, -0.5, -1/6]
    i == 28 && return SA[ 0.5, -0.5,  1/6]
    i == 29 && return SA[ 0.5,  0.5, -1/6]
    i == 30 && return SA[ 0.5,  0.5,  1/6]
    i == 31 && return SA[-0.5,  0.5, -1/6]
    i == 32 && return SA[-0.5,  0.5,  1/6]
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractHexahedron"))
end

function shapefunc_N(::Type{<:AbstractHexahedron}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    x, y, z = λ .* 2
    # Node dofs
    i == 1  && return 1/64 * (1 - x) * (1 - y) * (1 - z) * (9 * (x^2 + y^2 + z^2) - 19)
    i == 2  && return 1/64 * (1 + x) * (1 - y) * (1 - z) * (9 * (x^2 + y^2 + z^2) - 19)
    i == 3  && return 1/64 * (1 + x) * (1 + y) * (1 - z) * (9 * (x^2 + y^2 + z^2) - 19)
    i == 4  && return 1/64 * (1 - x) * (1 + y) * (1 - z) * (9 * (x^2 + y^2 + z^2) - 19)
    i == 5  && return 1/64 * (1 - x) * (1 - y) * (1 + z) * (9 * (x^2 + y^2 + z^2) - 19)
    i == 6  && return 1/64 * (1 + x) * (1 - y) * (1 + z) * (9 * (x^2 + y^2 + z^2) - 19)
    i == 7  && return 1/64 * (1 + x) * (1 + y) * (1 + z) * (9 * (x^2 + y^2 + z^2) - 19)
    i == 8  && return 1/64 * (1 - x) * (1 + y) * (1 + z) * (9 * (x^2 + y^2 + z^2) - 19)
    # Edge dofs
    i == 9  && return 9/64 * (1 - x^2) * (1 - 3 * x) * (1 - y) * (1 - z)
    i == 10 && return 9/64 * (1 - x^2) * (1 + 3 * x) * (1 - y) * (1 - z)
    i == 11 && return 9/64 * (1 - y^2) * (1 - 3 * y) * (1 + x) * (1 - z)
    i == 12 && return 9/64 * (1 - y^2) * (1 + 3 * y) * (1 + x) * (1 - z)
    i == 13 && return 9/64 * (1 - x^2) * (1 + 3 * x) * (1 + y) * (1 - z)
    i == 14 && return 9/64 * (1 - x^2) * (1 - 3 * x) * (1 + y) * (1 - z)
    i == 15 && return 9/64 * (1 - y^2) * (1 + 3 * y) * (1 - x) * (1 - z)
    i == 16 && return 9/64 * (1 - y^2) * (1 - 3 * y) * (1 - x) * (1 - z)

    i == 17 && return 9/64 * (1 - x^2) * (1 - 3 * x) * (1 - y) * (1 + z)
    i == 18 && return 9/64 * (1 - x^2) * (1 + 3 * x) * (1 - y) * (1 + z)
    i == 19 && return 9/64 * (1 - y^2) * (1 - 3 * y) * (1 + x) * (1 + z)
    i == 20 && return 9/64 * (1 - y^2) * (1 + 3 * y) * (1 + x) * (1 + z)
    i == 21 && return 9/64 * (1 - x^2) * (1 + 3 * x) * (1 + y) * (1 + z)
    i == 22 && return 9/64 * (1 - x^2) * (1 - 3 * x) * (1 + y) * (1 + z)
    i == 23 && return 9/64 * (1 - y^2) * (1 + 3 * y) * (1 - x) * (1 + z)
    i == 24 && return 9/64 * (1 - y^2) * (1 - 3 * y) * (1 - x) * (1 + z)

    i == 25 && return 9/64 * (1 - z^2) * (1 - 3 * z) * (1 - x) * (1 - y)
    i == 26 && return 9/64 * (1 - z^2) * (1 + 3 * z) * (1 - x) * (1 - y)
    i == 27 && return 9/64 * (1 - z^2) * (1 - 3 * z) * (1 + x) * (1 - y)
    i == 28 && return 9/64 * (1 - z^2) * (1 + 3 * z) * (1 + x) * (1 - y)
    i == 29 && return 9/64 * (1 - z^2) * (1 - 3 * z) * (1 + x) * (1 + y)
    i == 30 && return 9/64 * (1 - z^2) * (1 + 3 * z) * (1 + x) * (1 + y)
    i == 31 && return 9/64 * (1 - z^2) * (1 - 3 * z) * (1 - x) * (1 + y)
    i == 32 && return 9/64 * (1 - z^2) * (1 + 3 * z) * (1 - x) * (1 + y)
    throw(ArgumentError("No dof $i for Lagrang3{3} shape functions in AbstractHexahedron"))
end


function shapefunc_∇N(::Type{<:AbstractHexahedron}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    x, y, z = λ .* 2
    # Node dofs
    i == 1  && return SA[-1/32 * (1 - y) * (1 - z) * (9 * (3 * x^2 - 2 * x + y^2 + z^2) - 19),
                         -1/32 * (1 - x) * (1 - z) * (9 * (3 * y^2 - 2 * y + x^2 + z^2) - 19),
                         -1/32 * (1 - x) * (1 - y) * (9 * (3 * z^2 - 2 * z + x^2 + y^2) - 19)]
    i == 2  && return SA[ 1/32 * (1 - y) * (1 - z) * (9 * (3 * x^2 + 2 * x + y^2 + z^2) - 19),
                         -1/32 * (1 + x) * (1 - z) * (9 * (3 * y^2 - 2 * y + x^2 + z^2) - 19),
                         -1/32 * (1 + x) * (1 - y) * (9 * (3 * z^2 - 2 * z + x^2 + y^2) - 19)]
    i == 3  && return SA[ 1/32 * (1 + y) * (1 - z) * (9 * (3 * x^2 + 2 * x + y^2 + z^2) - 19),
                          1/32 * (1 + x) * (1 - z) * (9 * (3 * y^2 + 2 * y + x^2 + z^2) - 19),
                         -1/32 * (1 + x) * (1 + y) * (9 * (3 * z^2 - 2 * z + x^2 + y^2) - 19)]
    i == 4  && return SA[-1/32 * (1 + y) * (1 - z) * (9 * (3 * x^2 - 2 * x + y^2 + z^2) - 19),
                          1/32 * (1 - x) * (1 - z) * (9 * (3 * y^2 + 2 * y + x^2 + z^2) - 19),
                         -1/32 * (1 - x) * (1 + y) * (9 * (3 * z^2 - 2 * z + x^2 + y^2) - 19)]
    i == 5  && return SA[-1/32 * (1 - y) * (1 + z) * (9 * (3 * x^2 - 2 * x + y^2 + z^2) - 19),
                         -1/32 * (1 - x) * (1 + z) * (9 * (3 * y^2 - 2 * y + x^2 + z^2) - 19),
                          1/32 * (1 - x) * (1 - y) * (9 * (3 * z^2 + 2 * z + x^2 + y^2) - 19)]
    i == 6  && return SA[ 1/32 * (1 - y) * (1 + z) * (9 * (3 * x^2 + 2 * x + y^2 + z^2) - 19),
                         -1/32 * (1 + x) * (1 + z) * (9 * (3 * y^2 - 2 * y + x^2 + z^2) - 19),
                          1/32 * (1 + x) * (1 - y) * (9 * (3 * z^2 + 2 * z + x^2 + y^2) - 19)]
    i == 7  && return SA[ 1/32 * (1 + y) * (1 + z) * (9 * (3 * x^2 + 2 * x + y^2 + z^2) - 19),
                          1/32 * (1 + x) * (1 + z) * (9 * (3 * y^2 + 2 * y + x^2 + z^2) - 19),
                          1/32 * (1 + x) * (1 + y) * (9 * (3 * z^2 + 2 * z + x^2 + y^2) - 19)]
    i == 8  && return SA[-1/32 * (1 + y) * (1 + z) * (9 * (3 * x^2 - 2 * x + y^2 + z^2) - 19),
                          1/32 * (1 - x) * (1 + z) * (9 * (3 * y^2 + 2 * y + x^2 + z^2) - 19),
                          1/32 * (1 - x) * (1 + y) * (9 * (3 * z^2 + 2 * z + x^2 + y^2) - 19)]
    # Edge dofs
    i == 9  && return SA[ 9/32 * (1 - y) * (1 - z) * (9 * x^2 - 2 * x - 3),
                         -9/32 * (1 - x^2) * (1 - 3 * x) * (1 - z),
                         -9/32 * (1 - x^2) * (1 - 3 * x) * (1 - y)]
    i == 10 && return SA[-9/32 * (1 - y) * (1 - z) * (9 * x^2 + 2 * x - 3),
                         -9/32 * (1 - x^2) * (1 + 3 * x) * (1 - z),
                         -9/32 * (1 - x^2) * (1 + 3 * x) * (1 - y)]
    i == 11 && return SA[ 9/32 * (1 - y^2) * (1 - 3 * y) * (1 - z),
                          9/32 * (1 + x) * (1 - z) * (9 * y^2 - 2 * y - 3),
                         -9/32 * (1 - y^2) * (1 - 3 * y) * (1 + x)]
    i == 12 && return SA[ 9/32 * (1 - y^2) * (1 + 3 * y) * (1 - z),
                         -9/32 * (1 + x) * (1 - z) * (9 * y^2 + 2 * y - 3),
                         -9/32 * (1 - y^2) * (1 + 3 * y) * (1 + x)]
    i == 13 && return SA[-9/32 * (1 + y) * (1 - z) * (9 * x^2 + 2 * x - 3),
                          9/32 * (1 - x^2) * (1 + 3 * x) * (1 - z),
                         -9/32 * (1 - x^2) * (1 + 3 * x) * (1 + y)]
    i == 14 && return SA[ 9/32 * (1 + y) * (1 - z) * (9 * x^2 - 2 * x - 3),
                          9/32 * (1 - x^2) * (1 - 3 * x) * (1 - z),
                         -9/32 * (1 - x^2) * (1 - 3 * x) * (1 + y)]
    i == 15 && return SA[-9/32 * (1 - y^2) * (1 + 3 * y) * (1 - z),
                         -9/32 * (1 - x) * (1 - z) * (9 * y^2 + 2 * y - 3),
                         -9/32 * (1 - y^2) * (1 + 3 * y) * (1 - x)]
    i == 16 && return SA[-9/32 * (1 - y^2) * (1 - 3 * y) * (1 - z),
                          9/32 * (1 - x) * (1 - z) * (9 * y^2 - 2 * y - 3),
                         -9/32 * (1 - y^2) * (1 - 3 * y) * (1 - x)]

    i == 17 && return SA[ 9/32 * (1 - y) * (1 + z) * (9 * x^2 - 2 * x - 3),
                         -9/32 * (1 - x^2) * (1 - 3 * x) * (1 + z),
                          9/32 * (1 - x^2) * (1 - 3 * x) * (1 - y)]
    i == 18 && return SA[-9/32 * (1 - y) * (1 + z) * (9 * x^2 + 2 * x - 3),
                         -9/32 * (1 - x^2) * (1 + 3 * x) * (1 + z),
                          9/32 * (1 - x^2) * (1 + 3 * x) * (1 - y)]
    i == 19 && return SA[ 9/32 * (1 - y^2) * (1 - 3 * y) * (1 + z),
                          9/32 * (1 + x) * (1 + z) * (9 * y^2 - 2 * y - 3),
                          9/32 * (1 - y^2) * (1 - 3 * y) * (1 + x)]
    i == 20 && return SA[ 9/32 * (1 - y^2) * (1 + 3 * y) * (1 + z),
                         -9/32 * (1 + x) * (1 + z) * (9 * y^2 + 2 * y - 3),
                          9/32 * (1 - y^2) * (1 + 3 * y) * (1 + x)]
    i == 21 && return SA[-9/32 * (1 + y) * (1 + z) * (9 * x^2 + 2 * x - 3),
                          9/32 * (1 - x^2) * (1 + 3 * x) * (1 + z),
                          9/32 * (1 - x^2) * (1 + 3 * x) * (1 + y)]
    i == 22 && return SA[ 9/32 * (1 + y) * (1 + z) * (9 * x^2 - 2 * x - 3),
                          9/32 * (1 - x^2) * (1 - 3 * x) * (1 + z),
                          9/32 * (1 - x^2) * (1 - 3 * x) * (1 + y)]
    i == 23 && return SA[-9/32 * (1 - y^2) * (1 + 3 * y) * (1 + z),
                         -9/32 * (1 - x) * (1 + z) * (9 * y^2 + 2 * y - 3),
                          9/32 * (1 - y^2) * (1 + 3 * y) * (1 - x)]
    i == 24 && return SA[-9/32 * (1 - y^2) * (1 - 3 * y) * (1 + z),
                          9/32 * (1 - x) * (1 + z) * (9 * y^2 - 2 * y - 3),
                          9/32 * (1 - y^2) * (1 - 3 * y) * (1 - x)]

    i == 25 && return SA[-9/32 * (1 - z^2) * (1 - 3 * z) * (1 - y),
                         -9/32 * (1 - z^2) * (1 - 3 * z) * (1 - x),
                          9/32 * (1 - x) * (1 - y) * (9 * z^2 - 2 * z - 3)]
    i == 26 && return SA[-9/32 * (1 - z^2) * (1 + 3 * z) * (1 - y),
                         -9/32 * (1 - z^2) * (1 + 3 * z) * (1 - x),
                         -9/32 * (1 - x) * (1 - y) * (9 * z^2 + 2 * z - 3)]
    i == 27 && return SA[ 9/32 * (1 - z^2) * (1 - 3 * z) * (1 - y),
                         -9/32 * (1 - z^2) * (1 - 3 * z) * (1 + x),
                          9/32 * (1 + x) * (1 - y) * (9 * z^2 - 2 * z - 3)]
    i == 28 && return SA[ 9/32 * (1 - z^2) * (1 + 3 * z) * (1 - y),
                         -9/32 * (1 - z^2) * (1 + 3 * z) * (1 + x),
                         -9/32 * (1 + x) * (1 - y) * (9 * z^2 + 2 * z - 3)]
    i == 29 && return SA[ 9/32 * (1 - z^2) * (1 - 3 * z) * (1 + y),
                          9/32 * (1 - z^2) * (1 - 3 * z) * (1 + x),
                          9/32 * (1 + x) * (1 + y) * (9 * z^2 - 2 * z - 3)]
    i == 30 && return SA[ 9/32 * (1 - z^2) * (1 + 3 * z) * (1 + y),
                          9/32 * (1 - z^2) * (1 + 3 * z) * (1 + x),
                         -9/32 * (1 + x) * (1 + y) * (9 * z^2 + 2 * z - 3)]
    i == 31 && return SA[-9/32 * (1 - z^2) * (1 - 3 * z) * (1 + y),
                          9/32 * (1 - z^2) * (1 - 3 * z) * (1 - x),
                          9/32 * (1 - x) * (1 + y) * (9 * z^2 - 2 * z - 3)]
    i == 32 && return SA[-9/32 * (1 - z^2) * (1 + 3 * z) * (1 + y),
                          9/32 * (1 - z^2) * (1 + 3 * z) * (1 - x),
                         -9/32 * (1 - x) * (1 + y) * (9 * z^2 + 2 * z - 3)]
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractHexahedron"))
end

function edgedofs(::Type{<:AbstractHexahedron}, ::Lagrange{3}, i::Integer)
    i == 1  && return SA[1, 2,  9, 10]
    i == 2  && return SA[2, 3, 11, 12]
    i == 3  && return SA[3, 4, 13, 14]
    i == 4  && return SA[4, 1, 15, 16]
    i == 5  && return SA[5, 6, 17, 18]
    i == 6  && return SA[6, 7, 19, 20]
    i == 7  && return SA[7, 8, 21, 22]
    i == 8  && return SA[8, 5, 23, 24]
    i == 9  && return SA[1, 5, 25, 26]
    i == 10 && return SA[2, 6, 27, 28]
    i == 11 && return SA[3, 7, 29, 30]
    i == 12 && return SA[4, 8, 31, 32]
    throw(ArgumentError("No local edge $i in AbstractHexahedron"))
end

function facedofs(::Type{<:AbstractHexahedron}, ::Lagrange{3}, i::Integer)
    i == 1 && return SA[1, 2, 6, 5, 9, 10, 27, 28, 18, 17, 26, 25]
    i == 2 && return SA[2, 3, 7, 6, 11, 12, 29, 30, 20, 19, 28, 27]
    i == 3 && return SA[3, 4, 8, 7, 13, 14, 31, 32, 22, 21, 30, 29]
    i == 4 && return SA[4, 1, 5, 8, 15, 16, 25, 26, 24, 23, 32, 31]
    i == 5 && return SA[1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10,  9]
    i == 6 && return SA[5, 6, 7, 8, 17, 18, 19, 20, 21, 22, 23, 24]
    throw(ArgumentError("No local face $i in AbstractHexahedron"))
end

##############################
# Wedge
##############################

center_refcoordinates(::Type{<:AbstractWedge}) = SA[1/3, 1/3, 0.0]

function edge2cellcoordinates(::Type{<:AbstractWedge}, i::Integer, ξ::AbstractVector)
    i ==  1 && return SA[ ξ[1] + 0.5, 0.0, -0.5]
    i ==  2 && return SA[-ξ[1] + 0.5, ξ[1] + 0.5, -0.5]
    i ==  3 && return SA[0.0, -ξ[1] + 0.5, -0.5]
    i ==  4 && return SA[ ξ[1] + 0.5, 0.0, 0.5]
    i ==  5 && return SA[-ξ[1] + 0.5, ξ[1] + 0.5, 0.5]
    i ==  6 && return SA[0.0, -ξ[1] + 0.5, 0.5]
    i ==  7 && return SA[0.0, 0.0, ξ[1]]
    i ==  8 && return SA[1.0, 0.0, ξ[1]]
    i ==  9 && return SA[0.0, 1.0, ξ[1]]
    throw(ArgumentError("No local edge $i in shape type AbstractWedge"))
end

function face2cellcoordinates(::Type{<:AbstractWedge}, i::Integer, ξ::AbstractVector)
    i == 1 && return SA[ ξ[1] + 0.5, 0.0, ξ[2]]
    i == 2 && return SA[-ξ[1] + 0.5, ξ[1] + 0.5, ξ[2]]
    i == 3 && return SA[0.0, -ξ[1] + 0.5, ξ[2]]
    i == 4 && return SA[ξ[2], ξ[1], -0.5]
    i == 5 && return SA[ξ[1], ξ[2],  0.5]
    throw(ArgumentError("No local edge $i in shape type AbstractWedge"))
end

function vertexdofs(::Type{<:AbstractWedge}, ::Lagrange, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    i == 5 && return 5
    i == 6 && return 6
    throw(ArgumentError("No local vertex $i in AbstractWedge"))
end

function refcoordinates_grid(::Type{<:AbstractWedge}, n)
    # Length = sum(1:n) * n
    return (SA[ξ, η, ζ]  for ζ in LinRange(-0.5, 0.5, n)
        for (i, η) in enumerate(LinRange(0.0, 1.0, n))
            for ξ in LinRange(0.0, 1.0 - η, n - i + 1))
end

# --------------------------- Polynomial Order 1 ---------------------------
#=
Linear Wedge defined by 6 points

                  ○ (0,1,0.5)
                 /│6 ⋱
                / │     ⋱
               /  │        ⋱
    (0,0,0.5) ○───────────────○ (1,0,0.5)
              │4  │           │5
              │   │           │
              │   │           │
              │ 3 ○ (0,1,-0.5)│
              │  /   ⋱        │
              │ /       ⋱     │
              │/           ⋱  │
            1 ○───────────────○ 2
        (0,0,-0.5)        (1,0,-0.5)
=#

ndofs(       ::Type{<:AbstractWedge}, ::Lagrange{1}) = 6
ndofs_vertex(::Type{<:AbstractWedge}, ::Lagrange{1}) = (1, 1, 1, 1, 1, 1)
ndofs_edge(  ::Type{<:AbstractWedge}, ::Lagrange{1}) = (0, 0, 0, 0, 0, 0, 0, 0, 0)
ndofs_face(  ::Type{<:AbstractWedge}, ::Lagrange{1}) = (0, 0, 0, 0, 0)
ndofs_cell(  ::Type{<:AbstractWedge}, ::Lagrange{1}) = (0,)

function dof_refcoordinates(::Type{<:AbstractWedge}, ::Lagrange{1}, i::Integer)
    # Node dofs
    i == 1 && return SA[0.0, 0.0, -0.5]
    i == 2 && return SA[1.0, 0.0, -0.5]
    i == 3 && return SA[0.0, 1.0, -0.5]
    i == 4 && return SA[0.0, 0.0,  0.5]
    i == 5 && return SA[1.0, 0.0,  0.5]
    i == 6 && return SA[0.0, 1.0,  0.5]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractWedge"))
end

function shapefunc_N(::Type{<:AbstractWedge}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    x, y, z = λ
    # Node dofs
    i == 1 && return (1 - x - y) * (0.5 - z)
    i == 2 && return x * (0.5 - z)
    i == 3 && return y * (0.5 - z)
    i == 4 && return (1 - x - y) * (0.5 + z)
    i == 5 && return x * (0.5 + z)
    i == 6 && return y * (0.5 + z)
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractWedge"))
end

function shapefunc_∇N(::Type{<:AbstractWedge}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    x, y, z = λ
    # Node dofs
    i == 1  && return SA[-(0.5 - z),
                         -(0.5 - z),
                         -(1 - x - y)]
    i == 2  && return SA[ (0.5 - z),
                          0.0,
                         -x]
    i == 3  && return SA[ 0.0,
                          (0.5 - z),
                         -y]
    i == 4  && return SA[-(0.5 + z),
                         -(0.5 + z),
                          (1 - x - y)]
    i == 5  && return SA[ (0.5 + z),
                          0.0,
                          x]
    i == 6  && return SA[ 0.0,
                          (0.5 + z),
                          y]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractWedge"))
end

function edgedofs(::Type{<:AbstractWedge}, ::Lagrange{1}, i::Integer)
    i == 1  && return SA[1, 2]
    i == 2  && return SA[2, 3]
    i == 3  && return SA[3, 1]
    i == 4  && return SA[4, 5]
    i == 5  && return SA[5, 6]
    i == 6  && return SA[6, 4]
    i == 7  && return SA[1, 4]
    i == 8  && return SA[2, 5]
    i == 9  && return SA[3, 6]
    throw(ArgumentError("No local edge $i in AbstractWedge"))
end

# Type instability!
function facedofs(::Type{<:AbstractWedge}, ::Lagrange{1}, i::Integer)
    i == 1 && return SA[1, 2, 5, 4]
    i == 2 && return SA[2, 3, 6, 5]
    i == 3 && return SA[3, 1, 4, 6]
    i == 4 && return SA[1, 3, 2]
    i == 5 && return SA[4, 5, 6]
    throw(ArgumentError("No local face $i in AbstractWedge"))
end

# --------------------------- Polynomial Order 2 ---------------------------

ndofs(       ::Type{<:AbstractWedge}, ::Lagrange{2}) = 15
ndofs_vertex(::Type{<:AbstractWedge}, ::Lagrange{2}) = (1, 1, 1, 1, 1, 1)
ndofs_edge(  ::Type{<:AbstractWedge}, ::Lagrange{2}) = (1, 1, 1, 1, 1, 1, 1, 1, 1)
ndofs_face(  ::Type{<:AbstractWedge}, ::Lagrange{2}) = (0, 0, 0, 0, 0)
ndofs_cell(  ::Type{<:AbstractWedge}, ::Lagrange{2}) = (0,)

function dof_refcoordinates(::Type{<:AbstractWedge}, ::Lagrange{2}, i::Integer)
    # Node dofs
    i == 1  && return SA[0.0, 0.0, -0.5]
    i == 2  && return SA[1.0, 0.0, -0.5]
    i == 3  && return SA[0.0, 1.0, -0.5]
    i == 4  && return SA[0.0, 0.0,  0.5]
    i == 5  && return SA[1.0, 0.0,  0.5]
    i == 6  && return SA[0.0, 1.0,  0.5]
    # Edge dofs
    i == 7  && return SA[0.5, 0.0, -0.5]
    i == 8  && return SA[0.5, 0.5, -0.5]
    i == 9  && return SA[0.0, 0.5, -0.5]
    i == 10 && return SA[0.5, 0.0,  0.5]
    i == 11 && return SA[0.5, 0.5,  0.5]
    i == 12 && return SA[0.0, 0.5,  0.5]
    i == 13 && return SA[0.0, 0.0,  0.0]
    i == 14 && return SA[1.0, 0.0,  0.0]
    i == 15 && return SA[0.0, 1.0,  0.0]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractWedge"))
end

function shapefunc_N(::Type{<:AbstractWedge}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    x, y, z = λ
    ξ = 1 - x - y
    # Node dofs
    i == 1  && return ξ * (ξ - 0.5) * (1 - 2 * z) - 0.5 * ξ * (1 - 4 * z^2)
    i == 2  && return x * (x - 0.5) * (1 - 2 * z) - 0.5 * x * (1 - 4 * z^2)
    i == 3  && return y * (y - 0.5) * (1 - 2 * z) - 0.5 * y * (1 - 4 * z^2)
    i == 4  && return ξ * (ξ - 0.5) * (1 + 2 * z) - 0.5 * ξ * (1 - 4 * z^2)
    i == 5  && return x * (x - 0.5) * (1 + 2 * z) - 0.5 * x * (1 - 4 * z^2)
    i == 6  && return y * (y - 0.5) * (1 + 2 * z) - 0.5 * y * (1 - 4 * z^2)
    # Edge dofs
    i == 7  && return ξ * x * (2.0 - 4 * z)
    i == 8  && return x * y * (2.0 - 4 * z)
    i == 9  && return y * ξ * (2.0 - 4 * z)
    i == 10 && return ξ * x * (2.0 + 4 * z)
    i == 11 && return x * y * (2.0 + 4 * z)
    i == 12 && return y * ξ * (2.0 + 4 * z)
    i == 13 && return ξ * (1 - 4 * z^2)
    i == 14 && return x * (1 - 4 * z^2)
    i == 15 && return y * (1 - 4 * z^2)
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractWedge"))
end

function shapefunc_∇N(::Type{<:AbstractWedge}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    x, y, z = λ
    ξ = 1 - x - y
    # Node dofs
    i == 1  && return SA[(-2 * ξ + 0.5) * (1 - 2 * z) + 0.5 * (1 - 4 * z^2),
                         (-2 * ξ + 0.5) * (1 - 2 * z) + 0.5 * (1 - 4 * z^2),
                          -2 * ξ * (ξ - 0.5) + 4 * ξ * z]
    i == 2  && return SA[( 2 * x - 0.5) * (1 - 2 * z) - 0.5 * (1 - 4 * z^2),
                          0.0,
                         -2 * x * (x - 0.5) + 4 * x * z]
    i == 3  && return SA[ 0.0,
                         ( 2 * y - 0.5) * (1 - 2 * z) - 0.5 * (1 - 4 * z^2),
                         -2 * y * (y - 0.5) + 4 * y * z]
    i == 4  && return SA[(-2 * ξ + 0.5) * (1 + 2 * z) + 0.5 * (1 - 4 * z^2),
                         (-2 * ξ + 0.5) * (1 + 2 * z) + 0.5 * (1 - 4 * z^2),
                         2 * ξ * (ξ - 0.5) + 4 * ξ * z]
    i == 5  && return SA[( 2 * x - 0.5) * (1 + 2 * z) - 0.5 * (1 - 4 * z^2),
                          0.0,
                         2 * x * (x - 0.5) + 4 * x * z]
    i == 6  && return SA[ 0.0,
                         ( 2 * y - 0.5) * (1 + 2 * z) - 0.5 * (1 - 4 * z^2),
                         2 * y * (y - 0.5) + 4 * y * z]
    # Edge dofs
    i == 7  && return SA[(ξ - x) * (2.0 - 4 * z),
                         -x * (2.0 - 4 * z),
                         -4.0 * ξ * x]
    i == 8  && return SA[ y * (2.0 - 4 * z),
                          x * (2.0 - 4 * z),
                         -4.0 * x * y]
    i == 9  && return SA[ -y * (2.0 - 4 * z),
                         (-y + ξ) * (2.0 - 4 * z),
                         -4.0 * y * ξ]
    i == 10 && return SA[(ξ - x) * (2.0 + 4 * z),
                         -x * (2.0 + 4 * z),
                          4.0 * ξ * x]
    i == 11 && return SA[ y * (2.0 + 4 * z),
                          x * (2.0 + 4 * z),
                          4.0 * x * y]
    i == 12 && return SA[ -y * (2.0 + 4 * z),
                         (-y + ξ) * (2.0 + 4 * z),
                          4.0 * y * ξ]
    i == 13 && return SA[-(1.0 - 4 * z^2),
                         -(1.0 - 4 * z^2),
                         -8.0 * z * (1 - x - y)]
    i == 14 && return SA[ 1.0 - 4 * z^2,
                          0.0,
                         -8.0 * x * z]
    i == 15 && return SA[ 0.0,
                          1.0 - 4 * z^2,
                         -8.0 * y * z]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractWedge"))
end

function edgedofs(::Type{<:AbstractWedge}, ::Lagrange{2}, i::Integer)
    i == 1  && return SA[1, 2, 7]
    i == 2  && return SA[2, 3, 8]
    i == 3  && return SA[3, 1, 9]
    i == 4  && return SA[4, 5, 10]
    i == 5  && return SA[5, 6, 11]
    i == 6  && return SA[6, 4, 12]
    i == 7  && return SA[1, 4, 13]
    i == 8  && return SA[2, 5, 14]
    i == 9  && return SA[3, 6, 15]
    throw(ArgumentError("No local edge $i in AbstractWedge"))
end

# Type instability!
function facedofs(::Type{<:AbstractWedge}, ::Lagrange{2}, i::Integer)
    i == 1 && return SA[1, 2, 5, 4, 7, 14, 10, 13]
    i == 2 && return SA[2, 3, 6, 5, 8, 15, 11, 14]
    i == 3 && return SA[3, 1, 4, 6, 9, 13, 12, 15]
    i == 4 && return SA[1, 3, 2, 9, 8, 7]
    i == 5 && return SA[4, 5, 6, 10, 11, 12]
    throw(ArgumentError("No local face $i in AbstractWedge"))
end
