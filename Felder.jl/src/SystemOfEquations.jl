export SystemOfEquations
export NonSymmetric
export Symmetric_L
export Symmetric_U
export preassemble_square
export preassemble_rectangular
export preassemble_blocks
export update_solution_fields!

abstract type SymmetryType end

struct NonSymmetric <: SymmetryType end
struct Symmetric_U <: SymmetryType end
struct Symmetric_L <: SymmetryType end

struct SystemOfEquations{TA, T, S<:SymmetryType}
    jacobian::TA
    residual::Vector{T}
    symmetry::S
    Δu::Vector{T}

    csr_cols::Vector{Vector{Int}}

    # fieldnames::Vector{Symbol}
    fieldoffsets::Vector{Int}

    dirichletindices::Vector{Int}
end

function Base.show(io::Core.IO, ::MIME"text/plain", s::SystemOfEquations)
    println(io, typeof(s), " J Δu = r with Jacobian J:")
    display(s.jacobian)
    # println(io, "  Field names:                ", s.fieldnames)
    println(io, "  Field offsets:              ", s.fieldoffsets)
    println(io, "  J (min, max):               ", extrema(s.jacobian))
    println(io, "  r (min, max):               ", extrema(s.residual))
    println(io, "  Δu (min, max):              ", extrema(s.Δu))
end

function Base.show(io::Core.IO, ::MIME"text/plain", s::SystemOfEquations{<:Any, <:Complex})
    println(io, typeof(s), " J Δu = r with Jacobian J:")
    display(s.jacobian)
    # println(io, "  Field names:                ", s.fieldnames)
    println(io, "  Field offsets:              ", s.fieldoffsets)
end

symmetrytype(::SystemOfEquations{TA, T, S}) where {TA, T, S} = S
maxdof(dofs) = maximum(maximum, dofs[i] for i in eachindex(dofs) if isassigned(dofs, i))

"""
"""
preassemble_soe(T::Type{<:Number}, args...) = preassemble_soe(SparseMatrixCSC, T, args...)
preassemble_soe(blockpattern::Matrix, args...) = preassemble_soe(SparseMatrixCSC, Float64, blockpattern, args...)
preassemble_soe(maindofs::Vector, args...) = preassemble_soe(SparseMatrixCSC, Float64, Matrix(I, length(maindofs), length(maindofs)), maindofs, args...)

function preassemble_soe(::Type{TM}, ::Type{T},
        blockpattern::Matrix, maindofs::Vector, ncomps::Vector=ones(Int, length(maindofs)),
        symmetry::SymmetryType=NonSymmetric()) where {TM<:AbstractMatrix, T<:Number}
    J = preassemble_blocks(TM, T, blockpattern, maindofs, ncomps, symmetry)
    r = zeros(T, size(J, 1))
    Δu = zeros(T, size(J, 1))
    csr_cols = get_csr_cols(J)
    fieldoffsets = cumsum([0; maxdof.(maindofs) .* ncomps])
    return SystemOfEquations{typeof(J), eltype(r), typeof(symmetry)}(J, r, symmetry, Δu, csr_cols, fieldoffsets, Int[])
end

"""
"""
function preassemble(pdes::Union{AbstractPDE, Tuple{Vararg{AbstractPDE}}}, bchandler::BCHandler, symmetry=NonSymmetric())
    soe = preassemble(pdes, symmetry)
    bcs = unique_dirichlet_bcs(bchandler)
    for (n, pde) in enumerate(pdes)
        field = pde.field
        fieldoffset = soe.fieldoffsets[n]
        indices = find_dirichlet_indices(field, bcs)
        for i in indices
            push!(soe.dirichletindices, fieldoffset + i)
        end
    end
    return soe
end

"""
"""
function preassemble(pdes::Union{AbstractPDE, Tuple{Vararg{AbstractPDE}}}, symmetry::SymmetryType=NonSymmetric())
    @assert all(ncomps(pde) == ncomps(pde.field) for pde in pdes) # Obsolete?

    T = promote_type((dtype(pde.field) for pde in pdes)...)
    blockpattern = Matrix(I, length(pdes), length(pdes)) # TODO: Off-diagonal for multiphysics
    dofs = [pde.field.dofs for pde in pdes]
    n_comps = [ncomps(pde.field) for pde in pdes]

    soe = preassemble_soe(SparseMatrixCSC, T, blockpattern, dofs, n_comps, symmetry)

    # TODO: This is a hack to get the fieldoffsets. Should be done in a better way.
    for (i, pde) in enumerate(pdes)
        pde.field.dofoffset[] = soe.fieldoffsets[i]
    end

    return preassemble_soe(SparseMatrixCSC, T, blockpattern, dofs, n_comps, symmetry)
end

#################################
# SparseArrays - SparseMatrixCSC
#################################

function preassemble(::Type{TA}, ::Type{T}, elementdofs; symmetry=NonSymmetric(), ncomps=1, fieldoffsets=[0]) where {TA<:SparseMatrixCSC, T<:Number}
    J = _preassemble_sparse(TA, T, symmetry, elementdofs, ncomps)
    r = zeros(T, size(J, 1))
    Δu = zeros(T, size(J, 1))
    csr_cols = get_csr_cols(J)
    return SystemOfEquations{typeof(J), eltype(r), typeof(symmetry)}(J, r, symmetry, Δu, csr_cols, fieldoffsets, Int[])
end

"""
"""
preassemble_square(T::Type{<:Number}, args...) = preassemble_square(SparseMatrixCSC, T, args...)
preassemble_square(dofs, args...) = preassemble_square(SparseMatrixCSC, Float64, dofs, args...)

function preassemble_square(::Type{<:SparseMatrixCSC}, T::Type{<:Number}, dofs, ncomps::Int=1, ::S=NonSymmetric()) where {S <: SymmetryType}
    N = maxdof(dofs) * ncomps
    _rows = [Int[] for _ in 1:N]
    for _dofs in dofs
        for j in _dofs
            for n in 1:ncomps
                j_n = (j - 1) * ncomps + n
                for i in _dofs
                    for m in 1:ncomps
                        i_m = (i - 1) * ncomps + m
                        if S <: NonSymmetric
                            unique_sorted_insert!(_rows[j_n], i_m)
                        elseif S <: Symmetric_L
                            if i_m >= j_n
                                unique_sorted_insert!(_rows[j_n], i_m)
                            end
                        elseif S <: Symmetric_U
                            if i_m <= j_n
                                unique_sorted_insert!(_rows[j_n], i_m)
                            end
                        else
                            error("Symmetry type $S not supported.")
                        end
                    end
                end
            end
        end
    end

    colptr = ones(Int, N + 1)
    rowval = Int[]
    for j in eachindex(_rows)
        colptr[j + 1] = colptr[j] + length(_rows[j])
        append!(rowval, _rows[j])
    end

    nzval = zeros(T, length(rowval))

    return SparseMatrixCSC{T, Int}(N, N, colptr, rowval, nzval)
end

"""
"""
preassemble_rectangular(T::Type{<:Number}, args...) = preassemble_rectangular(SparseMatrixCSC, T, args...)
preassemble_rectangular(idofs, jdofs, args...) = preassemble_rectangular(SparseMatrixCSC, Float64, idofs, jdofs, args...)

function preassemble_rectangular(::Type{<:SparseMatrixCSC}, T::Type{<:Number}, idofs, jdofs, i_ncomps::Int=1, j_ncomps::Int=1)
    @assert count(i -> true, idofs)   == count(i -> true, jdofs) "`idofs and `jdofs` must have the same length."
    M = maxdof(idofs) * i_ncomps
    N = maxdof(jdofs) * j_ncomps
    _rows = [Int[] for _ in 1:N]
    for (_idofs, _jdofs) in zip(idofs, jdofs)
        for j in _jdofs
            for n in 1:j_ncomps
                j_n = (j - 1) * j_ncomps + n
                for i in _idofs
                    for m in 1:i_ncomps
                        i_m = (i - 1) * i_ncomps + m
                        unique_sorted_insert!(_rows[j_n], i_m)
                    end
                end
            end
        end
    end

    colptr = ones(Int, N + 1)
    rowval = Int[]
    for j in eachindex(_rows)
        colptr[j + 1] = colptr[j] + length(_rows[j])
        append!(rowval, _rows[j])
    end

    nzval = zeros(T, length(rowval))

    return SparseMatrixCSC{T, Int}(M, N, colptr, rowval, nzval)
end

"""
"""
preassemble_blocks(T::Type{<:Number}, args...) = preassemble_blocks(SparseMatrixCSC, T, args...)
preassemble_blocks(blockpattern::Matrix, args...) = preassemble_blocks(SparseMatrixCSC, Float64, blockpattern, args...)
preassemble_blocks(maindofs::Vector, args...) = preassemble_blocks(SparseMatrixCSC, Float64, Matrix(I, length(maindofs), length(maindofs)), maindofs, args...)

function preassemble_blocks(::Type{<:SparseMatrixCSC}, T::Type{<:Number},
        blockpattern::Matrix, maindofs::Vector, ncomps::Vector=ones(Int, length(maindofs)),
        symmetry::SymmetryType=NonSymmetric())
    @assert size(blockpattern, 1) == size(blockpattern, 2)
    @assert size(blockpattern, 1) == length(maindofs)
    @assert length(maindofs) == length(ncomps)

    N = length(maindofs)
    maxdofs = [maxdof(dofs) for dofs in maindofs]

    if symmetry isa NonSymmetric
        Aij_flat = [blockpattern[i, j] == 1 ?
            preassemble_rectangular(SparseMatrixCSC, T, maindofs[i], maindofs[j], ncomps[i], ncomps[j]) :
            spzeros(T, maxdofs[i] * ncomps[i], maxdofs[j] * ncomps[j])
            for i in 1:N for j in 1:N]
        return sparse_hvcat(ntuple(i -> N, Val(N)), Aij_flat...)
    else
        @assert blockpattern == I "Symmetric preassembly only implemented for
            identity blockpattern (i.e. no off-diagonal entries)."
        Aii = [preassemble_square(SparseMatrixCSC, T, maindofs[i], ncomps[i], symmetry)
            for i in 1:N]
        return blockdiag(Aii...)
    end
end

function get_csr_cols(A::SparseMatrixCSC)
    N = size(A, 1)
    cols = [Int[] for i in 1:N]

    rows = rowvals(A)

    for j in axes(A, 2)
        for k in nzrange(A, j)
            row = rows[k]
            unique_sorted_insert!(cols[row], j)
        end
    end

    return cols
end

function zerocolumn!(A::SparseMatrixCSC, j::Integer)
    @boundscheck checkbounds(A, j, j)

    vals = nonzeros(A)

    for i in nzrange(A, j)
        @inbounds vals[i] = zero(eltype(A))
    end

    return A
end

function zerorow!(A::SparseMatrixCSC, i::Integer, connecteddofs=get_csr_cols(A)[i])
    @boundscheck checkbounds(A, i, i)
    @boundscheck checkbounds(A, connecteddofs, connecteddofs)

    rows = rowvals(A)
    vals = nonzeros(A)

    @inbounds for j in connecteddofs
        found = false
        for k in nzrange(A, j)
            row = rows[k]
            if row == i
                vals[k] = zero(eltype(A))
                found = true
                break
            elseif row > i
                break
            end
        end
        found || error("Entry for row $i in column $j not found in sparse matrix.")
    end
    return A
end

"""
TODO
"""
function enforce_dirichlet!(J, r, dirichletdofs, csr_cols=get_csr_cols(J))
    @boundscheck checkbounds(J, dirichletdofs, dirichletdofs)
    @boundscheck checkbounds(r, dirichletdofs)
    @boundscheck checkbounds(csr_cols, dirichletdofs)

    for dof in dirichletdofs
        zerocolumn!(J, dof)
        zerorow!(J, dof, csr_cols[dof])
        J[dof, dof] = one(eltype(J)) # = maximum(abs, J) for better condition number?
        r[dof] = zero(eltype(J)) # This should better be u - d?
    end

    return J, r
end

function enforce_dirichlet!(soe::SystemOfEquations)
    enforce_dirichlet!(soe.jacobian, soe.residual, soe.dirichletindices, soe.csr_cols)
    return soe
end

############################
# Misc
############################

function update_solution_fields!(pde::AbstractPDE, args...)
    update_solution_fields!(Tuple(pde), args...)
end

function update_solution_fields!(pdes::Tuple{Vararg{AbstractPDE, 1}}, soe::SystemOfEquations,
        damping=1.0)
    @assert length(soe.fieldoffsets) - 1 == 1
    @unpack Δu, fieldoffsets = soe
    pde1, = pdes
    pde1.field.u .-= damping .* reinterpret(eltype(pde1.field.u), Δu)
    return pdes
end

function update_solution_fields!(pdes::Tuple{Vararg{AbstractPDE, 2}}, soe::SystemOfEquations,
        damping=1.0)
    @assert length(soe.fieldoffsets) - 1 == 2
    Δu = soe.Δu
    offsets = soe.fieldoffsets

    pde1, pde2 = pdes
    pde1.field.u .-= damping .* reinterpret(eltype(pde1.field.u), @view Δu[(offsets[1] + 1):offsets[2]])
    pde2.field.u .-= damping .* reinterpret(eltype(pde2.field.u), @view Δu[(offsets[2] + 1):offsets[3]])
    return pdes
end

function update_solution_fields!(pdes::Tuple{Vararg{AbstractPDE, 3}}, soe::SystemOfEquations,
        damping=1.0)
    @assert length(soe.fieldoffsets) - 1 == 3
    Δu = soe.Δu
    offsets = soe.fieldoffsets

    pde1, pde2, pde3 = pdes
    pde1.field.u .-= damping .* reinterpret(eltype(pde1.field.u), @view Δu[(offsets[1] + 1):offsets[2]])
    pde2.field.u .-= damping .* reinterpret(eltype(pde2.field.u), @view Δu[(offsets[2] + 1):offsets[3]])
    pde3.field.u .-= damping .* reinterpret(eltype(pde3.field.u), @view Δu[(offsets[3] + 1):offsets[4]])
    return pdes
end

function update_solution_fields!(pdes::Tuple{Vararg{AbstractPDE, 4}}, soe::SystemOfEquations,
        damping=1.0)
    @assert length(soe.fieldoffsets) - 1 == 4
    Δu = soe.Δu
    offsets = soe.fieldoffsets

    pde1, pde2, pde3, pde4 = pdes
    pde1.field.u .-= damping .* reinterpret(eltype(pde1.field.u), @view Δu[(offsets[1] + 1):offsets[2]])
    pde2.field.u .-= damping .* reinterpret(eltype(pde2.field.u), @view Δu[(offsets[2] + 1):offsets[3]])
    pde3.field.u .-= damping .* reinterpret(eltype(pde3.field.u), @view Δu[(offsets[3] + 1):offsets[4]])
    pde4.field.u .-= damping .* reinterpret(eltype(pde4.field.u), @view Δu[(offsets[4] + 1):offsets[5]])
    return pdes
end

function update_solution_fields!(pdes::Tuple{Vararg{AbstractPDE, 5}}, soe::SystemOfEquations,
        damping=1.0)
    @assert length(soe.fieldoffsets) - 1 == 5
    Δu = soe.Δu
    offsets = soe.fieldoffsets

    pde1, pde2, pde3, pde4, pde5 = pdes
    pde1.field.u .-= damping .* reinterpret(eltype(pde1.field.u), @view Δu[(offsets[1] + 1):offsets[2]])
    pde2.field.u .-= damping .* reinterpret(eltype(pde2.field.u), @view Δu[(offsets[2] + 1):offsets[3]])
    pde3.field.u .-= damping .* reinterpret(eltype(pde3.field.u), @view Δu[(offsets[3] + 1):offsets[4]])
    pde4.field.u .-= damping .* reinterpret(eltype(pde4.field.u), @view Δu[(offsets[4] + 1):offsets[5]])
    pde5.field.u .-= damping .* reinterpret(eltype(pde5.field.u), @view Δu[(offsets[5] + 1):offsets[6]])
    return pdes
end
