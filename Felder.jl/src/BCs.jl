export AbstractBC
export AbstractDirichletBC
export AbstractIntegratedBC
export DirichletBC
export dirichlet!
export dirichletvars
export BCHandler
export applydirichlet!
export NeumannBC
export RobinBC
export bcvars

abstract type AbstractBC end

Base.length(::AbstractBC) = 1
Base.iterate(bc::AbstractBC) = (bc, nothing)
Base.iterate(::AbstractBC, ::Any) = nothing

###############################
# Dirichlet BC
###############################

abstract type AbstractDirichletBC <: AbstractBC end

struct DirichletBC{N, F1, F2} <: AbstractDirichletBC
    field::F2
    boundaries::Vector{Symbol}
    u::F1 # u(x, t)
    enabled::SVector{N, Bool}
end

"""
"""
function DirichletBC(field, boundarytags, u, enabled=trues(ncomps(field)))
    _btags = boundarytags isa Symbol ? [boundarytags] : boundarytags

    N = ncomps(field)

    @assert N == length(enabled)
    @assert N == length(evalargs(u, fill(NaN, ndims(field)), NaN))

    DirichletBC{N, typeof(u), typeof(field)}(field, _btags, u, enabled)
end

ncomps(::DirichletBC{N}) where {N} = N

function Base.show(io::Core.IO, ::MIME"text/plain", bc::DirichletBC)
    println(io, "DirichletBC{$(ncomps(bc)), ...} for:")
    println(io, "  Field:            $(typeof(bc.field).name.wrapper){$(eltype(bc.field)), ...}")
    println(io, "  Boundaries:       $(bc.boundaries)")
    println(io, "  Prescribed value: $(bc.u)")
    print(io,   "  Enabled:          $(bc.enabled)")
end

function dirichletvars(bc::DirichletBC, fields, proxy, q)
    # This would require calling set_lagrange_evaluation!, then q
    # would be the dof location on the boundary
    return NamedTuple()
end

@inline function dirichlet!(u, bc::DirichletBC, vars, proxy, q)
    x = proxy.x[q]
    t = bc.field.t[]

    u .= evalargs(bc.u, x, t)
end

# -----------------------------------------------------------------

function find_dirichlet_indices(field::AbstractField, bcs)
    indices = Set{Int}()
    for bc in bcs
        if bc isa DirichletBC && bc.field === field
            _find_dirichlet_indices!(indices, bc)
        end
    end
    return sort!(collect(indices))
end

function _find_dirichlet_indices!(indices, bc::AbstractDirichletBC)
    field = bc.field
    shapefunc = field.shapefunc
    mesh = field.mesh
    dofs = field.dofs
    n_comps = ncomps(field)

    for facetindex in eachfacetindex(mesh, bc.boundaries)
        shape = mesh.elements[facetindex.elementid]
        shapedofs = dofs[facetindex.elementid]
        localfacetdofs = facetdofs(typeof(shape), shapefunc, facetindex.ilocal)

        for i in localfacetdofs
            dof = shapedofs[i]
            for j in 1:n_comps
                if bc.enabled[j]
                    push!(indices, (dof - 1) * n_comps + j)
                end
            end
        end
    end

    return indices
end

###############################
# Integrated BCs
###############################

abstract type AbstractIntegratedBC <: AbstractBC end

"""
    NeumannBC(field, boundarytags, g)

Adds a boundary contribution of the weak form

    ∫ a(x, t) Nᵢ dΓ

to the system of equations of `field`, where Γ is the boundary
defined by `boundarytags` and w are the weak form test functions.

`a` can be a Number or (Static) vector, or a function of the
coordinates `x` and time `t`.

Note: This simply adds the boundary contribution as described above.
The mathematical and physical constraint that is implied for a field
depends on the weak form of the PDE.

TODO: Specialize and rename this for PDE types to make sure the
condition/constraint is understood.
"""
struct NeumannBC{F1, F2} <: AbstractIntegratedBC
    field::F2
    boundaries::Vector{Symbol}
    a::F1 # g(x, t)
end

function NeumannBC(field, boundarytags, a)
    _btags = boundarytags isa Symbol ? [boundarytags] : boundarytags

    NeumannBC{typeof(a), typeof(field)}(field, _btags, a)
end

function Base.show(io::Core.IO, ::MIME"text/plain", bc::NeumannBC)
    println(io, "NeumannBC{...} contribution ∫ a(x, t) Nᵢ dΓ for:")
    println(io, "  Field:            $(typeof(bc.field).name.wrapper){$(eltype(bc.field)), ...}")
    println(io, "  Boundaries:       $(bc.boundaries)")
    print(io,   "  a:                $(bc.a)")
end

function bcvars(bc::NeumannBC, fields, proxy, q)
    x = proxy.x[q]
    t = bc.field.t[]
    return (
        a = evalargs(bc.a, x, t),
        x = x,
        t = t,
        N = proxy.N[q],
    )
end

@inline function jacobian!(J, ::NeumannBC, vars, proxy, i, j)
end

@inline function residual!(r, ::NeumannBC, vars, proxy, i)
    @unpack a, N = vars
    r[1] = a * N[i]
end

# -----------------------------------------------------------------

"""
    RobinBC(field, boundarytags, g, a)

Adds a boundary contribution of the weak form

    ∫ (a(x, t) + b(x, t) u) Nᵢ dΓ

to the system of equations of `field`, where Γ is the boundary
defined by `boundarytags` and w are the weak form test functions.

`a` and `b` can be a Number or (Static) vector, or a function of
the coordinates `x` and time `t`.

Note: This simply adds the boundary contribution as described above.
The mathematical and physical constraint that is implied for a field
depends on the weak form of the PDE.

TODO: Specialize and rename this for PDE types to make sure the
condition/constraint is understood.
"""
struct RobinBC{F1, F2, F3} <: AbstractIntegratedBC
    field::F3
    boundaries::Vector{Symbol}
    a::F1 # a(x, t)
    b::F2 # b(x, t)
end

function RobinBC(field, boundarytags, a, b)
    _btags = boundarytags isa Symbol ? [boundarytags] : boundarytags

    RobinBC{typeof(a), typeof(b), typeof(field)}(field, _btags, a, b)
end

function Base.show(io::Core.IO, ::MIME"text/plain", bc::RobinBC)
    println(io, "RobinBC{...} contribution ∫ (a + b u) Nᵢ dΓ for:")
    println(io, "  Field:            $(typeof(bc.field).name.wrapper){$(eltype(bc.field)), ...}")
    println(io, "  Boundaries:       $(bc.boundaries)")
    println(io, "  a:                $(bc.a)")
    print(io,   "  b:                $(bc.b)")
end

function bcvars(bc::RobinBC, fields, proxy, q)
    x = proxy.x[q]
    t = bc.field.t[]
    return (
        u = interpolate(bc.field, q),
        a = evalargs(bc.a, x, t),
        b = evalargs(bc.b, x, t),
        x = x,
        t = t,
        N = proxy.N[q],
    )
end

@inline function jacobian!(J, ::RobinBC, vars, proxy, i, j)
    @unpack b, N = vars
    J[1, 1] = b * N[j] * N[i]
end

@inline function residual!(r, ::RobinBC, vars, proxy, i)
    @unpack a, b, u, N = vars
    r[1] = (a + b * u) * N[i]
end

###############################
# BC Handler
###############################

struct BCHandler{M}
    dirichlet::Dict{Symbol, Tuple{Vararg{AbstractDirichletBC}}}
    integrated::Dict{Symbol, Tuple{Vararg{AbstractIntegratedBC}}}
    mesh::M
end

function BCHandler()
    dirichlet = Dict{Symbol, Tuple{Vararg{AbstractDirichletBC}}}()
    integrated = Dict{Symbol, Tuple{Vararg{AbstractIntegratedBC}}}()
    BCHandler{Nothing}(dirichlet, integrated, nothing)
end

BCHandler(::Tuple{}) = BCHandler()

function BCHandler(mesh::AbstractMesh)
    dirichlet = Dict{Symbol, Tuple{Vararg{AbstractDirichletBC}}}()
    integrated = Dict{Symbol, Tuple{Vararg{AbstractIntegratedBC}}}()
    BCHandler{typeof(mesh)}(dirichlet, integrated, mesh)
end

function BCHandler(bcs::Union{AbstractBC, Tuple{Vararg{AbstractBC}}})
    _dirichlet = Dict{Symbol, Vector{AbstractDirichletBC}}()
    _integrated = Dict{Symbol, Vector{AbstractIntegratedBC}}()

    mesh = first(bc.field.mesh for bc in bcs)

    for bc in bcs
        @assert bc.field.mesh === mesh "All BCs passed to a BCHandler must be defined on the same mesh"
        if bc isa AbstractDirichletBC
            for tag in bc.boundaries
                if !haskey(_dirichlet, tag)
                    _dirichlet[tag] = AbstractDirichletBC[]
                end
                push!(_dirichlet[tag], bc)
            end
        elseif bc isa AbstractIntegratedBC
            for tag in bc.boundaries
                if !haskey(_integrated, tag)
                    _integrated[tag] = AbstractIntegratedBC[]
                end
                push!(_integrated[tag], bc)
            end
        else
            error("$(typeof(bc)) not supported in BCHandler")
        end
    end

    dirichlet = Dict(key => tuple(value...) for (key, value) in _dirichlet)
    integrated = Dict(key => tuple(value...) for (key, value) in _integrated)

    BCHandler{typeof(mesh)}(dirichlet, integrated, mesh)
end

function Base.show(io::Core.IO, ::MIME"text/plain", bch::BCHandler)
    print(io, "BCHandler handling:")
    for (tag, bcs) in bch.dirichlet
        print(io, "\n  $(length(bcs)) Dirichlet BC(s) on :$tag")
    end
    for (tag, bcs) in bch.integrated
        print(io, "\n  $(length(bcs)) integrated BC(s) on :$tag")
    end
end

function applydirichlet!(bchandler::BCHandler, fieldhandler; kwargs...)
    for (tag, bcs) in bchandler.dirichlet
        applydirichlet!(bchandler.mesh, tag, bcs, fieldhandler; kwargs...)
    end

    return
end

function eachdirichletbc(bchandler::BCHandler)
    # This loop is not type stable
    return ((eachfacetindex(bchandler.mesh, tag), bcs) for (tag, bcs) in bchandler.dirichlet)
end

function eachintegratedbc(bchandler::BCHandler)
    # This loop is not type stable
    return ((eachfacetindex(bchandler.mesh, tag), bcs) for (tag, bcs) in bchandler.integrated)
end

function unique_dirichlet_bcs(bchandler::BCHandler)
    # This loop is not type stable
    return tuple(unique(bc for (tag, bcs) in bchandler.dirichlet for bc in bcs)...)
end

function unique_integrated_bcs(bchandler::BCHandler)
    # This loop is not type stable
    return tuple(unique(bc for (tag, bcs) in bchandler.integrated for bc in bcs)...)
end

###############################
# Dirichlet Application
###############################

function preallocate(bcs::Union{AbstractBC, Tuple{Vararg{AbstractBC}}})
    T = promote_type((dtype(bc.field) for bc in bcs)...)
    n = maximum(ncomps(bc.field) for bc in bcs)
    u = zeros(T, n)
    return u
end

function applydirichlet!(mesh::AbstractMesh, boundarytags, args...; kwargs...)
    facetindices = eachfacetindex(mesh, boundarytags)

    applydirichlet!(facetindices, args...; kwargs...)
end

function applydirichlet!(facetindices, bcs, fieldhandler)
    fields = fieldhandler.fields
    ud = preallocate(bcs) # TODO: Multithreading?

    for facetindex in facetindices
        setproxy_dof_evaluation!(fields, facetindex)
        applydirichlet_element!(ud, bcs, fields)
    end
end

function applydirichlet_element!(ud, bcs::AbstractDirichletBC, args...; kwargs...)
    applydirichlet_element!(ud, (bcs,), args...; kwargs...)
end

function applydirichlet_element!(ud, bcs::Tuple{Vararg{AbstractDirichletBC, 1}}, fields)
    bc1, = bcs

    _applydirichlet_element!(ud, bc1, fields)
end

function applydirichlet_element!(ud, bcs::Tuple{Vararg{AbstractDirichletBC, 2}}, fields)
    bc1, bc2 = bcs

    _applydirichlet_element!(ud, bc1, fields)
    _applydirichlet_element!(ud, bc2, fields)
end

function applydirichlet_element!(ud, bcs::Tuple{Vararg{AbstractDirichletBC, 3}}, fields)
    bc1, bc2, bc3 = bcs

    _applydirichlet_element!(ud, bc1, fields)
    _applydirichlet_element!(ud, bc2, fields)
    _applydirichlet_element!(ud, bc3, fields)
end

function applydirichlet_element!(ud, bcs::Tuple{Vararg{AbstractDirichletBC, 4}}, fields)
    bc1, bc2, bc3, bc4 = bcs

    _applydirichlet_element!(ud, bc1, fields)
    _applydirichlet_element!(ud, bc2, fields)
    _applydirichlet_element!(ud, bc3, fields)
    _applydirichlet_element!(ud, bc4, fields)
end

function applydirichlet_element!(ud, bcs::Tuple{Vararg{AbstractDirichletBC, 5}}, fields)
    bc1, bc2, bc3, bc4, bc5 = bcs

    _applydirichlet_element!(ud, bc1, fields)
    _applydirichlet_element!(ud, bc2, fields)
    _applydirichlet_element!(ud, bc3, fields)
    _applydirichlet_element!(ud, bc4, fields)
    _applydirichlet_element!(ud, bc5, fields)
end


function _applydirichlet_element!(ud, bc::AbstractDirichletBC, fields, mat=UndefinedMaterial())
    u = bc.field.u
    proxy = bc.field.proxies[threadid()]
    dofs = proxy.dofs
    localdofs = proxy.localdofs

    for q in eachindex(proxy.localdofs)
        u_dof = u[dofs[localdofs[q]]]
        vars = dirichletvars(bc, fields, proxy, q) # TODO: time, material

        dirichlet!(ud, bc, vars, proxy, q)

        u[dofs[localdofs[q]]] = _convert_enabled(u_dof, ud, bc.enabled)
    end

    return ud
end

function _convert_enabled(u_dof::T, ud, enabled) where {T <: Number}
    @assert length(u_dof) == length(enabled)

    return T(enabled[1] ? ud[1] : u_dof)
end

function _convert_enabled(u_dof::T, ud, enabled) where {T <: AbstractVector}
    @assert length(u_dof) == length(enabled)

    return T(enabled[i] ? ud[i] : u_dof[i] for i in eachindex(u_dof))
end

