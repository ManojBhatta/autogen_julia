export preallocate
export assemble_element!
export assemble_element
export preassemble
export assemble!
export SingleThreaded
export MultiThreaded

######################################
# Sparse Matrix Assignment (Assembly)
######################################

"""
TODO
Equivalent to

    A[dofs, dofs] += Ae

with `size(Ae) == (length(dofs), length(dofs))`
"""
function _assemble!(A::SparseMatrixCSC, dofs, Ae::AbstractMatrix,
                    dofperm=sortperm(dofs, alg=InsertionSort),
                    n_comps=1, offset=0, symmetry=NonSymmetric())
    @assert length(dofs) == length(dofperm)

    n_dofs = length(dofs)
    @boundscheck @assert min(size(A)...) >= maximum(dofs) * n_comps
    @boundscheck @assert min(size(Ae)...) >= n_dofs * n_comps

    rows = rowvals(A)
    vals = nonzeros(A)

    @inbounds for col_counter in eachindex(dofperm)
        jp = dofperm[col_counter]
        j = (jp - 1) * n_comps
        J = (dofs[jp] - 1) * n_comps + offset
        for n in 1:n_comps
            j_n = j + n
            J_n = J + n
            if symmetry isa Symmetric_L
                row_counter = col_counter
                m_min = n
            else
                row_counter = 1
                m_min = 1
            end
            ip = dofperm[row_counter]
            i = (ip - 1) * n_comps
            I_m = (dofs[ip] - 1) * n_comps + m_min + offset
            found = false
            k = A.colptr[J_n]
            k_max = A.colptr[J_n + 1]
            while k <= k_max
                if rows[k] == I_m
                    found = true
                    for m in m_min:n_comps
                        i_m = i + m
                        vals[k] += Ae[i_m, j_n]
                        k += 1
                    end
                    if row_counter < n_dofs
                        row_counter += 1
                        m_min = 1
                        ip = dofperm[row_counter]
                        i = (ip - 1) * n_comps
                        I_m = (dofs[ip] - 1) * n_comps + m_min + offset
                        found = false
                    else
                        break
                    end
                else
                    k += 1
                end
            end
            found || error("index [$I_m, $J_n] not found in preallocated sparse matrix.")
        end
    end
    return A
end

"""
TODO
Equivalent to

    A[idofs, jdofs] += Ae

with `size(Ae) == (length(idofs), length(jdofs))`
"""
function _assemble!(A::SparseMatrixCSC, idofs, jdofs, Ae::AbstractMatrix,
                    idofperm=sortperm(idofs, alg=InsertionSort), jdofperm=sortperm(jdofs, alg=InsertionSort),
                    i_ncomps=1, j_ncomps=1, ioffset=0, joffset=0)

    @assert length(idofs) == length(idofperm)
    @assert length(jdofs) == length(jdofperm)

    n_idofs = length(idofs)
    n_jdofs = length(jdofs)
    @boundscheck @assert size(A, 1) >= maximum(idofs) * i_ncomps
    @boundscheck @assert size(A, 2) >= maximum(jdofs) * j_ncomps
    @boundscheck @assert size(Ae, 1) >= n_idofs * i_ncomps
    @boundscheck @assert size(Ae, 2) >= n_jdofs * j_ncomps

    rows = rowvals(A)
    vals = nonzeros(A)

    @inbounds for col_counter in eachindex(jdofperm)
        jp = jdofperm[col_counter]
        j = (jp - 1) * j_ncomps
        J = (jdofs[jp] - 1) * j_ncomps + joffset
        for n in 1:j_ncomps
            j_n = j + n
            J_n = J + n
            row_counter = 1
            ip = idofperm[row_counter]
            i = (ip - 1) * i_ncomps
            I_m = (idofs[ip] - 1) * i_ncomps + 1 + ioffset
            found = false
            k = A.colptr[J_n]
            k_max = A.colptr[J_n + 1]
            while k <= k_max
                if rows[k] == I_m
                    found = true
                    for m in 1:i_ncomps
                        i_m = i + m
                        vals[k] += Ae[i_m, j_n]
                        k += 1
                    end
                    if row_counter < n_idofs
                        row_counter += 1
                        ip = idofperm[row_counter]
                        i = (ip - 1) * i_ncomps
                        I_m = (idofs[ip] - 1) * i_ncomps + 1 + ioffset
                        found = false
                    else
                        break
                    end
                else
                    k += 1
                end
            end
            found || error("index [$I_m, $J_n] not found in preallocated sparse matrix.")
        end
    end
    return A
end

"""
"""
function _assemble!(r::AbstractVector, dofs, re::AbstractVector, n_comps=1, offset=0)
    n_dofs = length(dofs)

    @boundscheck @assert length(r) >= maximum(n_dofs) * n_comps
    @boundscheck @assert length(re) >= n_dofs * n_comps

    @inbounds for k in 1:n_dofs
        for n in 1:n_comps
            I = (dofs[k] - 1) * n_comps + n + offset
            i = (k - 1) * n_comps + n
            r[I] += re[i]
        end
    end
    return r
end

"""
"""
function _assemble!(soe::SystemOfEquations, element_contribution, proxies, n_comps=1, offset=0)
    th_id = threadid()
    J, r = soe.jacobian, soe.residual
    Je, re = element_contribution.jacobian[th_id], element_contribution.residual[th_id]
    proxy = proxies[th_id]
    _assemble!(J, proxy.dofs, Je, proxy.dofperm, n_comps, offset, soe.symmetry)
    _assemble!(r, proxy.dofs, re, n_comps, offset)
    return soe
end

"""
"""
function _assemble!(soe::SystemOfEquations, element_contribution, field::AbstractField)
    proxies = field.proxies
    n_comps = ncomps(field)
    offset = field.dofoffset[]
    _assemble!(soe, element_contribution, proxies, n_comps, offset)
end

#############################
# Element Preallocation
#############################

struct ElementContribution{T, S<:SymmetryType}
    jacobian::Vector{Matrix{T}}
    residual::Vector{Vector{T}}
    symmetry::S
end

"""
    preallocate(shapetypes, pdes, symmetry=NonSymmetric())
    preallocate(shapes, pdes, symmetry=NonSymmetric())

Allocates the smallest possible dense element matrix and RHS-vector for `shapes`
with `shapefunc`.
TODO: Infer symmetry from pdes.
"""
function preallocate(shapetypes, pdes, symmetry=NonSymmetric())
    ncomps_per_field = (ncomps(pde.field) for pde in pdes)
    sf_per_field = (pde.field.shapefunc for pde in pdes)

    max_indices_per_field = (maximum(S -> ndofs(S, sf), shapetypes) * n_comps
        for (n_comps, sf) in zip(ncomps_per_field, sf_per_field))

    T = promote_type((dtype(pde.field) for pde in pdes)...)
    n = maximum(max_indices_per_field)

    J = [zeros(T, n, n) for _ in 1:nthreads()]
    r = [zeros(T, n) for _ in 1:nthreads()]

    return ElementContribution{T, typeof(symmetry)}(J, r, symmetry)
end

function preallocate(shapes::Vector{T}, args...; kwargs...) where {T <: AbstractShape}
    if isconcretetype(T)
        return preallocate((T,), args...; kwargs...)
    else
        unique_shapetypes = unique(typeof(shape) for shape in shapes)
        return preallocate(unique_shapetypes, args...; kwargs...)
    end
end

#############################
# Element Assembly
#############################

"""
    assemble_element(pde, fields, t [, material]; symmetry=NonSymmetric())

Assemble the Jacobian matrix and residual vector for the given `pde` and `fields`.
All proxies need to be set to the correct element and quadrature point before calling
using `set_element_integration!`.
"""
function assemble_element(pde, args...; symmetry=NonSymmetric(), kwargs...)
    shapes = pde.field.mesh.elements
    element_contr = preallocate(shapes, pde, symmetry)
    pde_contr = preallocate(pdes)
    assemble_element!(element_contr, pde_contr, args...; kwargs...)
end

"""
    assemble_element!(element_contribution, pde_contribution, pde, fields, t [, material])

Assemble the Jacobian matrix and residual vector the preallocated `element_contribution`
and `pde_contribution`.
TODO

# See also:
assemble_element
preallocate

# Examples:
"""
function assemble_element!(element_contribution::ElementContribution{T, NonSymmetric},
                           pde_contribution::PDEContribution,
                           pde, fields, mat=UndefinedMaterial()) where {T}
    th_id = threadid()
    Je, re = element_contribution.jacobian[th_id], element_contribution.residual[th_id]
    Jp, rp = pde_contribution.jacobian[th_id], pde_contribution.residual[th_id]

    proxy = pde.field.proxies[th_id]

    n_dofs = ndofs(proxy)
    n_comps = ncomps(pde)

    # TODO: Boundschecks?
    # @assert all(field.proxies[threadid()].current_mode[] == 0 for (_, field) in eachfield(fields))

    fill!(Je, 0)
    fill!(re, 0)
    # TODO: Fill Jp and rp with zeros?

    for q in eachindex(proxy.dΩ)
        dΩ = proxy.dΩ[q]
        vars = weakformvars(pde, fields, proxy, q)
        for j in 1:n_dofs
            j_m = (j - 1) * n_comps
            for i in 1:n_dofs
                i_m = (i - 1) * n_comps
                jacobian!(Jp, pde, vars, proxy, i, j)
                for n in 1:n_comps
                    for m in 1:n_comps
                        @inbounds Je[i_m + m, j_m + n] += Jp[m, n] * dΩ
                    end
                end
            end
            residual!(rp, pde, vars, proxy, j)
            for n in 1:n_comps
                @inbounds re[j_m + n] += rp[n] * dΩ
            end
        end
    end

    return element_contribution
end

"""
"""
function assemble_element_bc!(element_contribution::ElementContribution{T, NonSymmetric},
                              pde_contribution::PDEContribution,
                              bc::AbstractIntegratedBC,
                              fields, mat=UndefinedMaterial()) where {T}
    th_id = threadid()
    Je, re = element_contribution.jacobian[th_id], element_contribution.residual[th_id]
    Jp, rp = pde_contribution.jacobian[th_id], pde_contribution.residual[th_id]

    proxy = bc.field.proxies[th_id]

    n_dofs = ndofs(proxy)
    n_comps = ncomps(bc.field)

    # TODO: Boundschecks?
    # @assert all(field.proxies[threadid()].current_mode[] == 0 for (_, field) in eachfield(fields))

    fill!(Je, 0)
    fill!(re, 0)
    fill!(Jp, 0)
    fill!(rp, 0)

    for q in eachindex(proxy.dΩ)
        dΩ = proxy.dΩ[q]
        vars = bcvars(bc, fields, proxy, q)
        for j in 1:n_dofs
            j_m = (j - 1) * n_comps
            for i in 1:n_dofs
                i_m = (i - 1) * n_comps
                jacobian!(Jp, bc, vars, proxy, i, j)
                for n in 1:n_comps
                    for m in 1:n_comps
                        @inbounds Je[i_m + m, j_m + n] += Jp[m, n] * dΩ
                    end
                end
            end
            residual!(rp, bc, vars, proxy, j)
            for n in 1:n_comps
                @inbounds re[j_m + n] += rp[n] * dΩ
            end
        end
    end

    return element_contribution
end

###############################
# System of Equations Assembly
###############################

abstract type MultithreadingType end
struct SingleThreaded <: MultithreadingType end
struct MultiThreaded <: MultithreadingType end

"""
TODO
"""
function assemble!(soe::SystemOfEquations, pdes, fieldhandler::AbstractFieldHandler, bchandler::BCHandler,
                   material=UndefinedMaterial(), quadrature=defaultquadrature(fieldhandler),
                   threading=MultiThreaded())
    fields = fieldhandler.fields
    field, = fields
    mesh = field.mesh

    element_contr = preallocate(mesh.elements, pdes)
    pde_contr = preallocate(pdes)

    fill!(soe.jacobian, 0)
    fill!(soe.residual, 0)

    for elementindices in mesh.color2elements
        assemble!(soe, element_contr, pde_contr, elementindices, pdes, fieldhandler, material, quadrature, threading)
    end

    for (facetindices, bcs) in eachintegratedbc(bchandler)
        assemble_bcs!(soe, element_contr, pde_contr,facetindices, fieldhandler, bcs, material, quadrature)
    end

    enforce_dirichlet!(soe)

    return soe
end

"""
TODO
"""
function assemble!(soe::SystemOfEquations, element_contr, pde_contr, elementindices, pdes, fieldhandler,
                   material, quadrature, ::MultiThreaded)
    fields = fieldhandler.fields

    @assert npdes(pdes) == length(soe.fieldoffsets) - 1

    @threads :static for elementindex in elementindices
        setproxy_integration!(fields, elementindex, quadrature)
        assemble_element!(soe, element_contr, pde_contr, pdes, fields, material)
    end

    return soe
end

function assemble!(soe::SystemOfEquations, element_contr, pde_contr, elementindices, pdes, fieldhandler,
                   material, quadrature, ::SingleThreaded)
    fields = fieldhandler.fields

    @assert npdes(pdes) == length(soe.fieldoffsets) - 1

    for elementindex in elementindices
        setproxy_integration!(fields, elementindex, quadrature)
        assemble_element!(soe, element_contr, pde_contr, pdes, fields, material)
    end

    return soe
end

function assemble_bcs!(soe::SystemOfEquations, element_contr, pde_contr, facetindices, fieldhandler, bcs,
                       material, quadrature)
    fields = fieldhandler.fields

    for facetindex in facetindices
        setproxy_integration!(fields, facetindex, quadrature)
        assemble_element_bcs!(soe, element_contr, pde_contr, bcs, fields, material)
    end

    return soe
end

#################################
# Dispatch on PDEs
#################################
# TODO: generated functions

function assemble_element!(soe::SystemOfEquations, element_contr, pde_contr,
                           pde::AbstractPDE, args...; kwargs...)
    assemble_element!(soe, element_contr, pde_contr, (pde,), args...; kwargs...)
end

function assemble_element!(soe::SystemOfEquations, element_contr, pde_contr,
                           pdes::Tuple{Vararg{AbstractPDE, 1}}, fields, mat)
    pde1, = pdes

    assemble_element!(element_contr, pde_contr, pde1, fields, mat)
    _assemble!(soe, element_contr, pde1.field)

    return soe
end

function assemble_element!(soe::SystemOfEquations, element_contr, pde_contr,
                           pdes::Tuple{Vararg{AbstractPDE, 2}}, fields, mat)
    pde1, pde2 = pdes

    assemble_element!(element_contr, pde_contr, pde1, fields, mat)
    _assemble!(soe, element_contr, pde1.field)

    assemble_element!(element_contr, pde_contr, pde2, fields, mat)
    _assemble!(soe, element_contr, pde2.field)

    return soe
end

function assemble_element!(soe::SystemOfEquations, element_contr, pde_contr,
                           pdes::Tuple{Vararg{AbstractPDE, 3}}, fields, mat)
    pde1, pde2, pde3 = pdes

    assemble_element!(element_contr, pde_contr, pde1, fields, mat)
    _assemble!(soe, element_contr, pde1.field)

    assemble_element!(element_contr, pde_contr, pde2, fields, mat)
    _assemble!(soe, element_contr, pde2.field)

    assemble_element!(element_contr, pde_contr, pde3, fields, mat)
    _assemble!(soe, element_contr, pde3.field)

    return soe
end

function assemble_element!(soe::SystemOfEquations, element_contr, pde_contr,
                           pdes::Tuple{Vararg{AbstractPDE, 4}}, fields, mat)
    pde1, pde2, pde3, pde4 = pdes

    assemble_element!(element_contr, pde_contr, pde1, fields, mat)
    _assemble!(soe, element_contr, pde1.field)

    assemble_element!(element_contr, pde_contr, pde2, fields, mat)
    _assemble!(soe, element_contr, pde2.field)

    assemble_element!(element_contr, pde_contr, pde3, fields, mat)
    _assemble!(soe, element_contr, pde3.field)

    assemble_element!(element_contr, pde_contr, pde4, fields, mat)
    _assemble!(soe, element_contr, pde4.field)

    return soe
end

function assemble_element!(soe::SystemOfEquations, element_contr, pde_contr,
                           pdes::Tuple{Vararg{AbstractPDE, 5}}, fields, mat)
    pde1, pde2, pde3, pde4, pde5 = pdes

    assemble_element!(element_contr, pde_contr, pde1, fields, mat)
    _assemble!(soe, element_contr, pde1.field)

    assemble_element!(element_contr, pde_contr, pde2, fields, mat)
    _assemble!(soe, element_contr, pde2.field)

    assemble_element!(element_contr, pde_contr, pde3, fields, mat)
    _assemble!(soe, element_contr, pde3.field)

    assemble_element!(element_contr, pde_contr, pde4, fields, mat)
    _assemble!(soe, element_contr, pde4.field)

    assemble_element!(element_contr, pde_contr, pde5, fields, mat)
    _assemble!(soe, element_contr, pde5.field)

    return soe
end

function assemble_element_bcs!(soe::SystemOfEquations, element_contr, pde_contr,
                               bc::AbstractIntegratedBC, args...; kwargs...)
    assemble_element_bcs!(soe, element_contr, pde_contr, (bc,), args...; kwargs...)
end

function assemble_element_bcs!(soe::SystemOfEquations, element_contr, pde_contr,
                               bcs::Tuple{Vararg{AbstractIntegratedBC, 1}}, fields, mat)
    bc1, = bcs

    assemble_element_bc!(element_contr, pde_contr, bc1, fields, mat)
    _assemble!(soe, element_contr, bc1.field)

    return soe
end

function assemble_element_bcs!(soe::SystemOfEquations, element_contr, pde_contr,
                               bcs::Tuple{Vararg{AbstractIntegratedBC, 2}}, fields, mat)
    bc1, bc2 = bcs

    assemble_element_bc!(element_contr, pde_contr, bc1, fields, mat)
    _assemble!(soe, element_contr, bc1.field)

    assemble_element_bc!(element_contr, pde_contr, bc2, fields, mat)
    _assemble!(soe, element_contr, bc2.field)

    return soe
end

function assemble_element_bcs!(soe::SystemOfEquations, element_contr, pde_contr,
                               bcs::Tuple{Vararg{AbstractIntegratedBC, 3}}, fields, mat)
    bc1, bc2, bc3 = bcs

    assemble_element_bc!(element_contr, pde_contr, bc1, fields, mat)
    _assemble!(soe, element_contr, bc1.field)

    assemble_element_bc!(element_contr, pde_contr, bc2, fields, mat)
    _assemble!(soe, element_contr, bc2.field)

    assemble_element_bc!(element_contr, pde_contr, bc3, fields, mat)
    _assemble!(soe, element_contr, bc3.field)

    return soe
end

function assemble_element_bcs!(soe::SystemOfEquations, element_contr, pde_contr,
                               bcs::Tuple{Vararg{AbstractIntegratedBC, 4}}, fields, mat)
    bc1, bc2, bc3, bc4 = bcs

    assemble_element_bc!(element_contr, pde_contr, bc1, fields, mat)
    _assemble!(soe, element_contr, bc1.field)

    assemble_element_bc!(element_contr, pde_contr, bc2, fields, mat)
    _assemble!(soe, element_contr, bc2.field)

    assemble_element_bc!(element_contr, pde_contr, bc3, fields, mat)
    _assemble!(soe, element_contr, bc3.field)

    assemble_element_bc!(element_contr, pde_contr, bc4, fields, mat)
    _assemble!(soe, element_contr, bc4.field)

    return soe
end

function assemble_element_bcs!(soe::SystemOfEquations, element_contr, pde_contr,
                               bcs::Tuple{Vararg{AbstractIntegratedBC, 5}}, fields, mat)
    bc1, bc2, bc3, bc4, bc5 = bcs

    assemble_element_bc!(element_contr, pde_contr, bc1, fields, mat)
    _assemble!(soe, element_contr, bc1.field)

    assemble_element_bc!(element_contr, pde_contr, bc2, fields, mat)
    _assemble!(soe, element_contr, bc2.field)

    assemble_element_bc!(element_contr, pde_contr, bc3, fields, mat)
    _assemble!(soe, element_contr, bc3.field)

    assemble_element_bc!(element_contr, pde_contr, bc4, fields, mat)
    _assemble!(soe, element_contr, bc4.field)

    assemble_element_bc!(element_contr, pde_contr, bc5, fields, mat)
    _assemble!(soe, element_contr, bc5.field)

    return soe
end
