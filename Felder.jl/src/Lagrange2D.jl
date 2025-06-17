##############################
# Triangle
##############################

center_refcoordinates(::Type{<:AbstractTriangle}) = SA[1/3, 1/3]

function edge2facecoordinates(::Type{<:AbstractTriangle}, i::Integer, ξ::AbstractVector)
    i == 1 && return SA[ξ[1] + 0.5, 0.0]
    i == 2 && return SA[-ξ[1] + 0.5, ξ[1] + 0.5]
    i == 3 && return SA[0.0, -ξ[1] + 0.5]
    throw(ArgumentError("No local edge $i in shape type AbstractTriangle"))
end

function vertexdofs(::Type{<:AbstractTriangle}, ::Lagrange, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    throw(ArgumentError("No local vertex $i in AbstractTriangle"))
end

function refcoordinates_grid(::Type{<:AbstractTriangle}, n)
    # Length = sum(1:n)
    return (SA[ξ, η] for (i, η) in enumerate(LinRange(0.0, 1.0, n))
        for ξ in LinRange(0.0, 1.0 - η, n - i + 1))
end

# --------------------------- Polynomial Order 1 ---------------------------
#=
     (0,1)
     3 ○
       ┊ ⋱
       ┊   ⋱
       ┊     ⋱
       ┊       ⋱
     1 ○┄┄┄┄┄┄┄┄┄○ 2
     (0,0)       (1,0)
=#

ndofs(       ::Type{<:AbstractTriangle}, ::Lagrange{1}) = 3
ndofs_vertex(::Type{<:AbstractTriangle}, ::Lagrange{1}) = (1, 1, 1)
ndofs_edge(  ::Type{<:AbstractTriangle}, ::Lagrange{1}) = (0, 0, 0)
ndofs_face(  ::Type{<:AbstractTriangle}, ::Lagrange{1}) = (0,)

function dof_refcoordinates(::Type{<:AbstractTriangle}, ::Lagrange{1}, i::Integer)
    # Node dofs
    i == 1 && return SA[0.0, 0.0]
    i == 2 && return SA[1.0, 0.0]
    i == 3 && return SA[0.0, 1.0]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractTriangle"))
end

function shapefunc_N(::Type{<:AbstractTriangle}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    λ1, λ2 = λ
    λ3 = 1 - λ1 - λ2
    # Node dofs
    i == 1 && return λ3
    i == 2 && return λ1
    i == 3 && return λ2
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractTriangle"))
end

function shapefunc_∇N(::Type{<:AbstractTriangle}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    # Node dofs
    i == 1 && return SA[-1.0, -1.0]
    i == 2 && return SA[ 1.0,  0.0]
    i == 3 && return SA[ 0.0,  1.0]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractTriangle"))
end

function shapefunc_∇∇N(::Type{<:AbstractTriangle}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    # Node dofs
    i == 1 && return SA[0.0 0.0; 0.0 0.0]
    i == 2 && return SA[0.0 0.0; 0.0 0.0]
    i == 3 && return SA[0.0 0.0; 0.0 0.0]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractTriangle"))
end

function edgedofs(::Type{<:AbstractTriangle}, ::Lagrange{1}, i::Integer)
    i == 1 && return SA[1, 2]
    i == 2 && return SA[2, 3]
    i == 3 && return SA[3, 1]
    throw(ArgumentError("No local edge $i in AbstractTriangle"))
end

flip_facedofs_permutation(::Type{<:AbstractTriangle}, ::Lagrange{1}) = SVector{0, Int}()
circshift_facedofs_permutation(::Type{<:AbstractTriangle}, ::Lagrange{1}) = SVector{0, Int}()

# --------------------------- Polynomial Order 2 ---------------------------

ndofs(       ::Type{<:AbstractTriangle}, ::Lagrange{2}) = 6
ndofs_vertex(::Type{<:AbstractTriangle}, ::Lagrange{2}) = (1, 1, 1)
ndofs_edge(  ::Type{<:AbstractTriangle}, ::Lagrange{2}) = (1, 1, 1)
ndofs_face(  ::Type{<:AbstractTriangle}, ::Lagrange{2}) = (0,)

function dof_refcoordinates(::Type{<:AbstractTriangle}, ::Lagrange{2}, i::Integer)
    # Node dofs
    i == 1 && return SA[0.0, 0.0]
    i == 2 && return SA[1.0, 0.0]
    i == 3 && return SA[0.0, 1.0]
    # Edge dofs
    i == 4 && return SA[0.5, 0.0]
    i == 5 && return SA[0.5, 0.5]
    i == 6 && return SA[0.0, 0.5]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractTriangle"))
end

function shapefunc_N(::Type{<:AbstractTriangle}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    λ1, λ2 = λ
    λ3 = 1 - λ1 - λ2
    # Node dofs
    i == 1 && return λ3 * (2.0 * λ3 - 1)
    i == 2 && return λ1 * (2.0 * λ1 - 1)
    i == 3 && return λ2 * (2.0 * λ2 - 1)
    # Edge dofs
    i == 4 && return 4.0 * λ3 * λ1
    i == 5 && return 4.0 * λ1 * λ2
    i == 6 && return 4.0 * λ2 * λ3
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractTriangle"))
end

function shapefunc_∇N(::Type{<:AbstractTriangle}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    λ1, λ2 = λ
    λ3 = 1 - λ1 - λ2
    # Node dofs
    i == 1 && return SA[-4.0 * λ3 + 1,
                        -4.0 * λ3 + 1]
    i == 2 && return SA[ 4.0 * λ1 - 1,
                         0.0]
    i == 3 && return SA[ 0.0,
                         4.0 * λ2 - 1]
    # Edge dofs
    i == 4 && return SA[ 4.0 * (λ3 - λ1),
                        -4.0 * λ1]
    i == 5 && return SA[ 4.0 * λ2,
                         4.0 * λ1]
    i == 6 && return SA[-4.0 * λ2,
                         4.0 * (λ3 - λ2)]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractTriangle"))
end

function shapefunc_∇∇N(::Type{<:AbstractTriangle}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    λ1, λ2 = λ
    # Node dofs
    i == 1 && return SA[ 4.0  4.0;
                         4.0  4.0]
    i == 2 && return SA[ 4.0  0.0;
                         0.0  0.0]
    i == 3 && return SA[ 0.0  0.0;
                         0.0  4.0]
    # Edge dofs
    i == 4 && return SA[-8.0 -4.0;
                        -4.0  0.0]
    i == 5 && return SA[ 0.0  4.0;
                         4.0  0.0]
    i == 6 && return SA[ 0.0 -4.0;
                        -4.0 -8.0]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractTriangle"))
end

function edgedofs(::Type{<:AbstractTriangle}, ::Lagrange{2}, i::Integer)
    i == 1 && return SA[1, 2, 4]
    i == 2 && return SA[2, 3, 5]
    i == 3 && return SA[3, 1, 6]
    throw(ArgumentError("No local edge $i in AbstractTriangle"))
end

flip_facedofs_permutation(::Type{<:AbstractTriangle}, ::Lagrange{2}) = SVector{0, Int}()
circshift_facedofs_permutation(::Type{<:AbstractTriangle}, ::Lagrange{2}) = SVector{0, Int}()

# --------------------------- Polynomial Order 3 ---------------------------

ndofs(       ::Type{<:AbstractTriangle}, ::Lagrange{3}) = 10
ndofs_vertex(::Type{<:AbstractTriangle}, ::Lagrange{3}) = (1, 1, 1)
ndofs_edge(  ::Type{<:AbstractTriangle}, ::Lagrange{3}) = (2, 2, 2)
ndofs_face(  ::Type{<:AbstractTriangle}, ::Lagrange{3}) = (1,)

function dof_refcoordinates(::Type{<:AbstractTriangle}, ::Lagrange{3}, i::Integer)
    # Node dofs
    i == 1 && return SA[0.0, 0.0]
    i == 2 && return SA[1.0, 0.0]
    i == 3 && return SA[0.0, 1.0]
    # Edge dofs
    i == 4 && return SA[1/3, 0.0]
    i == 5 && return SA[2/3, 0.0]
    i == 6 && return SA[2/3, 1/3]
    i == 7 && return SA[1/3, 2/3]
    i == 8 && return SA[0.0, 2/3]
    i == 9 && return SA[0.0, 1/3]
    # Face dofs
    i == 10 && return SA[1/3, 1/3]
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractTriangle"))
end

function shapefunc_N(::Type{<:AbstractTriangle}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    λ1, λ2 = λ
    λ3 = 1 - λ1 - λ2
    # Node dofs
    i == 1 && return 0.5 * λ3 * (3 * λ3 - 1) * (3 * λ3 - 2)
    i == 2 && return 0.5 * λ1 * (3 * λ1 - 1) * (3 * λ1 - 2)
    i == 3 && return 0.5 * λ2 * (3 * λ2 - 1) * (3 * λ2 - 2)
    # Edge dofs
    i == 4 && return 4.5 * λ3 * λ1 * (3 * λ3 - 1)
    i == 5 && return 4.5 * λ3 * λ1 * (3 * λ1 - 1)
    i == 6 && return 4.5 * λ1 * λ2 * (3 * λ1 - 1)
    i == 7 && return 4.5 * λ1 * λ2 * (3 * λ2 - 1)
    i == 8 && return 4.5 * λ2 * λ3 * (3 * λ2 - 1)
    i == 9 && return 4.5 * λ2 * λ3 * (3 * λ3 - 1)
    # Face dofs
    i == 10 && return 27.0 * λ3 * λ1 * λ2
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractTriangle"))
end

function shapefunc_∇N(::Type{<:AbstractTriangle}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    λ1, λ2 = λ
    λ3 = 1 - λ1 - λ2
    # Node dofs
    i == 1  && return SA[-λ3 * (13.5 * λ3 - 9) - 1,
                         -λ3 * (13.5 * λ3 - 9) - 1]
    i == 2  && return SA[ λ1 * (13.5 * λ1 - 9) + 1,
                          0.0]
    i == 3  && return SA[ 0.0,
                          λ2 * (13.5 * λ2 - 9) + 1]
    # Edge dofs
    i == 4  && return SA[  4.5 * (λ3 * ( 3 * λ3 - 6 * λ1 - 1) + λ1),
                         -27.0 * λ1 * λ3 + 4.5 * λ1]
    i == 5  && return SA[  4.5 * (λ1 * (-3 * λ1 + 6 * λ3 + 1) - λ3),
                          -4.5 * λ1 * (3 * λ1 - 1)]
    i == 6  && return SA[ 27.0 * λ1 * λ2 - 4.5 * λ2,
                           4.5 * λ1 * (3 * λ1 - 1)]

    i == 7  && return SA[  4.5 * λ2 * (3 * λ2 - 1),
                          27.0 * λ2 * λ1 - 4.5 * λ1]
    i == 8  && return SA[ -4.5 * λ2 * (3 * λ2 - 1),
                           4.5 * (λ2 * (-3 * λ2 + 6 * λ3 + 1) - λ3)]
    i == 9  && return SA[-27.0 * λ2 * λ3 + 4.5 * λ2,
                           4.5 * (λ3 * ( 3 * λ3 - 6 * λ2 - 1) + λ2)]
    # Face dofs
    i == 10 && return SA[ 27.0 * λ2 * (λ3 - λ1),
                          27.0 * λ1 * (λ3 - λ2)]
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractTriangle"))
end

function edgedofs(::Type{<:AbstractTriangle}, ::Lagrange{3}, i::Integer)
    i == 1 && return SA[1, 2, 4, 5]
    i == 2 && return SA[2, 3, 6, 7]
    i == 3 && return SA[3, 1, 8, 9]
    throw(ArgumentError("No local edge $i in AbstractTriangle"))
end

flip_facedofs_permutation(::Type{<:AbstractTriangle}, ::Lagrange{3}) = SA[1]
circshift_facedofs_permutation(::Type{<:AbstractTriangle}, ::Lagrange{3}) =  SA[1]

##############################
# Quadrilateral
##############################

center_refcoordinates(::Type{<:AbstractQuadrilateral}) = SA[0.0, 0.0]

function edge2facecoordinates(::Type{<:AbstractQuadrilateral}, i::Integer, ξ::AbstractVector)
    i == 1 && return SA[ξ[1], -0.5]
    i == 2 && return SA[0.5, ξ[1]]
    i == 3 && return SA[-ξ[1], 0.5]
    i == 4 && return SA[-0.5, -ξ[1]]
    throw(ArgumentError("No local edge $i in shape type AbstractQuadrilateral"))
end

function vertexdofs(::Type{<:AbstractQuadrilateral}, ::Lagrange, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    i == 3 && return 3
    i == 4 && return 4
    throw(ArgumentError("No local vertex $i in AbstractQuadrilateral"))
end

function refcoordinates_grid(::Type{<:AbstractQuadrilateral}, n)
    # Length = n * n
    return (SA[ξ, η] for η in LinRange(-0.5, 0.5, n) for ξ in LinRange(-0.5, 0.5, n))
end

# --------------------------- Polynomial Order 1 ---------------------------
#=
    (-0.5,0.5)        (0.5,0.5)
        4 ○─────────────○ 3
          │             │
          │             │
          │             │
        1 ○─────────────○ 2
    (-0.5,-0.5)       (0.5,-0.5)
=#

ndofs(       ::Type{<:AbstractQuadrilateral}, ::Lagrange{1}) = 4
ndofs_vertex(::Type{<:AbstractQuadrilateral}, ::Lagrange{1}) = (1, 1, 1, 1)
ndofs_edge(  ::Type{<:AbstractQuadrilateral}, ::Lagrange{1}) = (0, 0, 0, 0)
ndofs_face(  ::Type{<:AbstractQuadrilateral}, ::Lagrange{1}) = (0,)

function dof_refcoordinates(::Type{<:AbstractQuadrilateral}, ::Lagrange{1}, i::Integer)
    # Node dofs
    i == 1 && return SA[-0.5, -0.5]
    i == 2 && return SA[ 0.5, -0.5]
    i == 3 && return SA[ 0.5,  0.5]
    i == 4 && return SA[-0.5,  0.5]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractQuadrilateral"))
end

function shapefunc_N(::Type{<:AbstractQuadrilateral}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    ξ, η = λ
    # Node dofs
    i == 1 && return (0.5 - ξ) * (0.5 - η)
    i == 2 && return (0.5 + ξ) * (0.5 - η)
    i == 3 && return (0.5 + ξ) * (0.5 + η)
    i == 4 && return (0.5 - ξ) * (0.5 + η)
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractQuadrilateral"))
end

function shapefunc_∇N(::Type{<:AbstractQuadrilateral}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    ξ, η = λ
    # Node dofs
    i == 1 && return SA[-0.5 + η,
                        -0.5 + ξ]
    i == 2 && return SA[ 0.5 - η,
                        -0.5 - ξ]
    i == 3 && return SA[ 0.5 + η,
                         0.5 + ξ]
    i == 4 && return SA[-0.5 - η,
                         0.5 - ξ]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractQuadrilateral"))
end

function shapefunc_∇∇N(::Type{<:AbstractQuadrilateral}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    # Node dofs
    i == 1 && return SA[ 0.0  1.0;
                         1.0  0.0]
    i == 2 && return SA[ 0.0 -1.0;
                        -1.0  0.0]
    i == 3 && return SA[ 0.0  1.0;
                         1.0  0.0]
    i == 4 && return SA[ 0.0 -1.0;
                        -1.0  0.0]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractQuadrilateral"))
end

function edgedofs(::Type{<:AbstractQuadrilateral}, ::Lagrange{1}, i::Integer)
    i == 1 && return SA[1, 2]
    i == 2 && return SA[2, 3]
    i == 3 && return SA[3, 4]
    i == 4 && return SA[4, 1]
    throw(ArgumentError("No local edge $i in AbstractQuadrilateral"))
end

flip_facedofs_permutation(::Type{<:AbstractQuadrilateral}, ::Lagrange{1}) = SVector{0, Int}()
circshift_facedofs_permutation(::Type{<:AbstractQuadrilateral}, ::Lagrange{1}) = SVector{0, Int}()

# --------------------------- Polynomial Order 2 ---------------------------
#=
  (-0.5,0.5)        (0.5,0.5)
       4 ○──────○──────○ 3
         │      7      │
         ○ 8         6 ○
         │      5      │
       1 ○──────○──────○ 2
  (-0.5,-0.5)       (0.5,-0.5)
=#

ndofs(       ::Type{<:AbstractQuadrilateral}, ::Lagrange{2}) = 8
ndofs_vertex(::Type{<:AbstractQuadrilateral}, ::Lagrange{2}) = (1, 1, 1, 1)
ndofs_edge(  ::Type{<:AbstractQuadrilateral}, ::Lagrange{2}) = (1, 1, 1, 1)
ndofs_face(  ::Type{<:AbstractQuadrilateral}, ::Lagrange{2}) = (0,)

function dof_refcoordinates(::Type{<:AbstractQuadrilateral}, ::Lagrange{2}, i::Integer)
    # Node dofs
    i == 1 && return SA[-0.5, -0.5]
    i == 2 && return SA[ 0.5, -0.5]
    i == 3 && return SA[ 0.5,  0.5]
    i == 4 && return SA[-0.5,  0.5]
    i == 5 && return SA[ 0.0, -0.5]
    i == 6 && return SA[ 0.5,  0.0]
    i == 7 && return SA[ 0.0,  0.5]
    i == 8 && return SA[-0.5,  0.0]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractQuadrilateral"))
end

function shapefunc_N(::Type{<:AbstractQuadrilateral}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    ξ, η = λ .* 2
    # Node dofs
    i == 1 && return 0.25 * (1 - ξ) * (1 - η) * (-1 - ξ - η)
    i == 2 && return 0.25 * (1 + ξ) * (1 - η) * (-1 + ξ - η)
    i == 3 && return 0.25 * (1 + ξ) * (1 + η) * (-1 + ξ + η)
    i == 4 && return 0.25 * (1 - ξ) * (1 + η) * (-1 - ξ + η)
    i == 5 && return 0.5  * (1 + ξ) * (1 - ξ) * (1 - η)
    i == 6 && return 0.5  * (1 + ξ) * (1 + η) * (1 - η)
    i == 7 && return 0.5  * (1 + ξ) * (1 - ξ) * (1 + η)
    i == 8 && return 0.5  * (1 - ξ) * (1 + η) * (1 - η)
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractQuadrilateral"))
end

function shapefunc_∇N(::Type{<:AbstractQuadrilateral}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    ξ, η = λ .* 2
    # Node dofs
    i == 1 && return SA[ξ * (1 - η) + 0.5 * η * (1 - η),
                        0.5 * ξ * (1 - ξ) + η * (1 - ξ)]
    i == 2 && return SA[ξ * ( 1 - η) + 0.5 * η * (-1 + η),
                        0.5 * ξ * (-1 - ξ) + η * ( 1 + ξ)]
    i == 3 && return SA[ ξ * (1 + η) + 0.5 * η * (1 + η),
                        0.5 * ξ * (1 + ξ) +  η * (1 + ξ)]
    i == 4 && return SA[ ξ * ( 1 + η) + 0.5 * η * (-1 - η),
                        0.5 * ξ * (-1 + ξ) +  η * ( 1 - ξ)]
    i == 5 && return SA[ 2 * ξ * (-1.0 + η),
                        -1.0 + ξ * ξ]
    i == 6 && return SA[ 1.0 - η * η,
                         2 * η * (-1.0 - ξ)]
    i == 7 && return SA[ 2 *  ξ * (-1.0 - η),
                         1.0 - ξ * ξ]
    i == 8 && return SA[-1.0 + η * η,
                         2 * η * (-1.0 + ξ)]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractQuadrilateral"))
end

function shapefunc_∇∇N(::Type{<:AbstractQuadrilateral}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    ξ, η = λ .* 2
    # Node dofs
    i == 1 && return 2.0 * SA[(1 - η)        (-ξ - η + 0.5);
                              (-ξ - η + 0.5) (1 - ξ)]
    i == 2 && return 2.0 * SA[(1 - η)        (-ξ + η - 0.5);
                              (-ξ + η - 0.5) (1 + ξ)]
    i == 3 && return 2.0 * SA[(1 + η)        ( ξ + η + 0.5);
                              ( ξ + η + 0.5) (1 + ξ)]
    i == 4 && return 2.0 * SA[(1 + η)        ( ξ - η - 0.5);
                              ( ξ - η - 0.5) (1 - ξ)]
    i == 5 && return 2.0 * SA[(-2 * (1 - η)) (2 * ξ);
                              (2 * ξ)        (0.0)]
    i == 6 && return 2.0 * SA[(0.0)          (-2 * η);
                              (-2 * η)       (-2 * (1 + ξ))]
    i == 7 && return 2.0 * SA[(-2 * (1 + η)) (-2 * ξ);
                              (-2 * ξ)       (0.0)]
    i == 8 && return 2.0 * SA[(0.0)          (2 * η);
                              (2 * η)        (-2 * (1 - ξ))]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractQuadrilateral"))
end

function edgedofs(::Type{<:AbstractQuadrilateral}, ::Lagrange{2}, i::Integer)
    i == 1 && return SA[1, 2, 5]
    i == 2 && return SA[2, 3, 6]
    i == 3 && return SA[3, 4, 7]
    i == 4 && return SA[4, 1, 8]
    throw(ArgumentError("No local edge $i in AbstractQuadrilateral"))
end

flip_facedofs_permutation(::Type{<:AbstractQuadrilateral}, ::Lagrange{2}) = SVector{0, Int}()
circshift_facedofs_permutation(::Type{<:AbstractQuadrilateral}, ::Lagrange{2}) = SVector{0, Int}()

# --------------------------- Polynomial Order 3 ---------------------------
#=
  (-0.5,0.5)        (0.5,0.5)
       4 ○───○─────○───○ 3
      11 ○   10    9   ○ 8
         │             │
      12 ○   5     6   ○ 7
       1 ○───○─────○───○ 2
  (-0.5,-0.5)       (0.5,-0.5)
=#

ndofs(       ::Type{<:AbstractQuadrilateral}, ::Lagrange{3}) = 12
ndofs_vertex(::Type{<:AbstractQuadrilateral}, ::Lagrange{3}) = (1, 1, 1, 1)
ndofs_edge(  ::Type{<:AbstractQuadrilateral}, ::Lagrange{3}) = (2, 2, 2, 2)
ndofs_face(  ::Type{<:AbstractQuadrilateral}, ::Lagrange{3}) = (0,)

function dof_refcoordinates(::Type{<:AbstractQuadrilateral}, ::Lagrange{3}, i::Integer)
    # Node dofs
    i ==  1 && return SA[-0.5, -0.5]
    i ==  2 && return SA[ 0.5, -0.5]
    i ==  3 && return SA[ 0.5,  0.5]
    i ==  4 && return SA[-0.5,  0.5]
    i ==  5 && return SA[-1/6, -0.5]
    i ==  6 && return SA[ 1/6, -0.5]
    i ==  7 && return SA[ 0.5, -1/6]
    i ==  8 && return SA[ 0.5,  1/6]
    i ==  9 && return SA[ 1/6,  0.5]
    i == 10 && return SA[-1/6,  0.5]
    i == 11 && return SA[-0.5,  1/6]
    i == 12 && return SA[-0.5, -1/6]
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractQuadrilateral"))
end

function shapefunc_N(::Type{<:AbstractQuadrilateral}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    ξ, η = λ .* 2
    # Node dofs
    i ==  1 && return 1/32 * (1 - ξ) * (1 - η) * (9 * (ξ^2 + η^2) - 10)
    i ==  2 && return 1/32 * (1 + ξ) * (1 - η) * (9 * (ξ^2 + η^2) - 10)
    i ==  3 && return 1/32 * (1 + ξ) * (1 + η) * (9 * (ξ^2 + η^2) - 10)
    i ==  4 && return 1/32 * (1 - ξ) * (1 + η) * (9 * (ξ^2 + η^2) - 10)
    i ==  5 && return 9/32 * (1 - ξ^2) * (1 - 3 * ξ) * (1 - η)
    i ==  6 && return 9/32 * (1 - ξ^2) * (1 + 3 * ξ) * (1 - η)
    i ==  7 && return 9/32 * (1 + ξ) * (1 - η^2) * (1 - 3 * η)
    i ==  8 && return 9/32 * (1 + ξ) * (1 - η^2) * (1 + 3 * η)
    i ==  9 && return 9/32 * (1 - ξ^2) * (1 + 3 * ξ) * (1 + η)
    i == 10 && return 9/32 * (1 - ξ^2) * (1 - 3 * ξ) * (1 + η)
    i == 11 && return 9/32 * (1 - ξ) * (1 - η^2) * (1 + 3 * η)
    i == 12 && return 9/32 * (1 - ξ) * (1 - η^2) * (1 - 3 * η)
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractQuadrilateral"))
end

function shapefunc_∇N(::Type{<:AbstractQuadrilateral}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    ξ, η = λ .* 2
    # Node dofs
    i ==  1 && return SA[-1/16 * (1 - η) * (9 * (3 * ξ^2 - 2 * ξ + η^2) - 10),
                         -1/16 * (1 - ξ) * (9 * (ξ^2 - 2 * η + 3 * η^2) - 10)]
    i ==  2 && return SA[ 1/16 * (1 - η) * (9 * (3 * ξ^2 + 2 * ξ + η^2) - 10),
                         -1/16 * (1 + ξ) * (9 * (ξ^2 - 2 * η + 3 * η^2) - 10)]
    i ==  3 && return SA[ 1/16 * (1 + η) * (9 * (3 * ξ^2 + 2 * ξ + η^2) - 10),
                          1/16 * (1 + ξ) * (9 * (ξ^2 + 2 * η + 3 * η^2) - 10)]
    i ==  4 && return SA[-1/16 * (1 + η) * (9 * (3 * ξ^2 - 2 * ξ + η^2) - 10),
                          1/16 * (1 - ξ) * (9 * (ξ^2 + 2 * η + 3 * η^2) - 10)]
    i ==  5 && return SA[ 9/16 * (1 - η) * (9 * ξ^2 - 2 * ξ - 3),
                         -9/16 * (1 - ξ^2) * (1 - 3 * ξ)]
    i ==  6 && return SA[-9/16 * (1 - η) * (9 * ξ^2 + 2 * ξ - 3),
                         -9/16 * (1 - ξ^2) * (1 + 3 * ξ)]
    i ==  7 && return SA[ 9/16 * (1 - η^2) * (1 - 3 * η),
                          9/16 * (1 + ξ) * (9 * η^2 - 2 * η - 3)]
    i ==  8 && return SA[ 9/16 * (1 - η^2) * (1 + 3 * η),
                         -9/16 * (1 + ξ) * (9 * η^2 + 2 * η - 3)]
    i ==  9 && return SA[-9/16 * (1 + η) * (9 * ξ^2 + 2 * ξ - 3),
                          9/16 * (1 - ξ^2) * (1 + 3 * ξ)]
    i == 10 && return SA[ 9/16 * (1 + η) * (9 * ξ^2 - 2 * ξ - 3),
                          9/16 * (1 - ξ^2) * (1 - 3 * ξ)]
    i == 11 && return SA[-9/16 * (1 - η^2) * (1 + 3 * η),
                        -9/16 * (1 - ξ) * (9 * η^2 + 2 * η - 3)]
    i == 12 && return SA[-9/16 * (1 - η^2) * (1 - 3 * η),
                          9/16 * (1 - ξ) * (9 * η^2 - 2 * η - 3)]
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractQuadrilateral"))
end

function edgedofs(::Type{<:AbstractQuadrilateral}, ::Lagrange{3}, i::Integer)
    i == 1 && return SA[1, 2, 5, 6]
    i == 2 && return SA[2, 3, 7, 8]
    i == 3 && return SA[3, 4, 9, 10]
    i == 4 && return SA[4, 1, 11, 12]
    throw(ArgumentError("No local edge $i in AbstractQuadrilateral"))
end

flip_facedofs_permutation(::Type{<:AbstractQuadrilateral}, ::Lagrange{3}) = SVector{0, Int}()
circshift_facedofs_permutation(::Type{<:AbstractQuadrilateral}, ::Lagrange{3}) = SVector{0, Int}()
