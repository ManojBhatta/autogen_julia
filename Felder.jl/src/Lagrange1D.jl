##############################
# Edge
##############################

center_refcoordinates(::Type{<:AbstractEdge}) = SA[0.0]

function vertex2edgecoordinates(::Type{<:AbstractEdge}, i::Integer, ξ::AbstractVector)
    i == 1 && return SA[-0.5]
    i == 2 && return SA[ 0.5]
    throw(ArgumentError("No local vertex $i in shape type AbstractEdge"))
end

function edge2edgecoordinates(::Type{<:AbstractEdge}, i::Integer, ξ::AbstractVector)
    i == 1 && return SA[ξ[1]]
    throw(ArgumentError("No local edge $i in shape type AbstractEdge"))
end

function vertexdofs(::Type{<:AbstractEdge}, ::Lagrange, i::Integer)
    i == 1 && return 1
    i == 2 && return 2
    throw(ArgumentError("No local vertex $i in AbstractEdge"))
end

function refcoordinates_grid(::Type{<:AbstractEdge}, n)
    # Length = n
    return (SA[ξ] for ξ in LinRange(-0.5, 0.5, n))
end

# --------------------------- Polynomial Order 1 ---------------------------
#=
      N1            N2
       ○────────────○
    -0.5            0.5
=#

ndofs(       ::Type{<:AbstractEdge}, ::Lagrange{1}) = 2
ndofs_vertex(::Type{<:AbstractEdge}, ::Lagrange{1}) = (1, 1)
ndofs_edge(  ::Type{<:AbstractEdge}, ::Lagrange{1}) = (0,)

function dof_refcoordinates(::Type{<:AbstractEdge}, ::Lagrange{1}, i::Integer)
    # Node dofs
    i == 1 && return SA[-0.5]
    i == 2 && return SA[ 0.5]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractEdge"))
end

function shapefunc_N(::Type{<:AbstractEdge}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    λ1, = λ
    # Node dofs
    i == 1 && return 0.5 - λ1
    i == 2 && return 0.5 + λ1
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractEdge"))
end

function shapefunc_∇N(::Type{<:AbstractEdge}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    # Node dofs
    i == 1 && return SA[-1.0]
    i == 2 && return SA[ 1.0]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractEdge"))
end

function shapefunc_∇∇N(::Type{<:AbstractEdge}, ::Lagrange{1}, i::Integer, λ::AbstractVector)
    # Node dofs
    i == 1 && return SA[0.0;;]
    i == 2 && return SA[0.0;;]
    throw(ArgumentError("No dof $i for Lagrange{1} shape functions in AbstractEdge"))
end

flip_edgedofs_permutation(::Type{<:AbstractEdge}, ::Lagrange{1}) = SVector{0, Int}()

# --------------------------- Polynomial Order 2 ---------------------------
#=
      N1     N3     N2
       ○──────○──────○
    -0.5      0      0.5
=#

ndofs(       ::Type{<:AbstractEdge}, ::Lagrange{2}) = 3
ndofs_vertex(::Type{<:AbstractEdge}, ::Lagrange{2}) = (1, 1)
ndofs_edge(  ::Type{<:AbstractEdge}, ::Lagrange{2}) = (1,)

function dof_refcoordinates(::Type{<:AbstractEdge}, ::Lagrange{2}, i::Integer)
    # Node dofs
    i == 1 && return SA[-0.5]
    i == 2 && return SA[ 0.5]
    # Edge dofs
    i == 3 && return SA[ 0.0]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractEdge"))
end

function shapefunc_N(::Type{<:AbstractEdge}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    λ1, = λ
    # Node dofs
    i == 1 && return -λ1 * (1 - 2.0 * λ1)
    i == 2 && return  λ1 * (1 + 2.0 * λ1)
    # Edge dofs
    i == 3 && return (1 + 2.0 * λ1) * (1 - 2.0 * λ1)
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractEdge"))
end

function shapefunc_∇N(::Type{<:AbstractEdge}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    λ1, = λ
    # Node dofs
    i == 1 && return SA[ 4.0 * λ1 - 1]
    i == 2 && return SA[ 4.0 * λ1 + 1]
    # Edge dofs
    i == 3 && return SA[-8.0 * λ1]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractEdge"))
end

function shapefunc_∇∇N(::Type{<:AbstractEdge}, ::Lagrange{2}, i::Integer, λ::AbstractVector)
    λ1, = λ
    # Node dofs
    i == 1 && return SA[ 4.0;;]
    i == 2 && return SA[ 4.0;;]
    # Edge dofs
    i == 3 && return SA[-8.0;;]
    throw(ArgumentError("No dof $i for Lagrange{2} shape functions in AbstractEdge"))
end

flip_edgedofs_permutation(::Type{<:AbstractEdge}, ::Lagrange{2}) = SA[1,]

# --------------------------- Polynomial Order 3 ---------------------------
#=
      N1   N3    N4   N2
       ○────○────○────○
    -0.5  -1/6  1/6   0.5
=#

ndofs(       ::Type{<:AbstractEdge}, ::Lagrange{3}) = 4
ndofs_vertex(::Type{<:AbstractEdge}, ::Lagrange{3}) = (1, 1)
ndofs_edge(  ::Type{<:AbstractEdge}, ::Lagrange{3}) = (2,)

function dof_refcoordinates(::Type{<:AbstractEdge}, ::Lagrange{3}, i::Integer)
    # Node dofs
    i == 1 && return SA[-0.5]
    i == 2 && return SA[ 0.5]
    # Edge dofs
    i == 3 && return SA[-1/6]
    i == 4 && return SA[ 1/6]
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractEdge"))
end

function shapefunc_N(::Type{<:AbstractEdge}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    λ1, = λ
    # Node dofs
    i == 1 && return -1/16 * (6 * λ1 + 1) * (6 * λ1 - 1) * (2 * λ1 - 1)
    i == 2 && return  1/16 * (6 * λ1 + 1) * (6 * λ1 - 1) * (2 * λ1 + 1)
    # Edge dofs
    i == 3 && return  9/16 * (2 * λ1 + 1) * (6 * λ1 - 1) * (2 * λ1 - 1)
    i == 4 && return -9/16 * (2 * λ1 + 1) * (6 * λ1 + 1) * (2 * λ1 - 1)
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractEdge"))
end

function shapefunc_∇N(::Type{<:AbstractEdge}, ::Lagrange{3}, i::Integer, λ::AbstractVector)
    λ1, = λ
    # Node dofs
    i == 1 && return SA[ 9/2 * λ1 * (-3 * λ1 + 1) + 1/8]
    i == 2 && return SA[ 9/2 * λ1 * ( 3 * λ1 + 1) - 1/8]
    # Edge dofs
    i == 3 && return SA[ 9/2 * λ1 * (9 * λ1 - 1) - 27/8]
    i == 4 && return SA[-9/2 * λ1 * (9 * λ1 + 1) + 27/8]
    throw(ArgumentError("No dof $i for Lagrange{3} shape functions in AbstractEdge"))
end

flip_edgedofs_permutation(::Type{<:AbstractEdge}, ::Lagrange{3}) = SA[2, 1]
