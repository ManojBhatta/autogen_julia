export ncomps
export evaluate
export AbstractField
export AbstractDomainField
export AbstractBoundaryField
export ScalarField
export ScalarElementalField
export VectorField
export VectorElementalField
export ScalarBoundaryField
export VectorBoundaryField
export eachfield
export eachshape
export eachshapeindex
export getdofmap
export FieldHandler
export FieldHandler!
export init!
export setproxy_integration!
export setproxy_dof_evaluation!
export setproxy_point_evaluation!
export initproxy_point_evaluation!
export remesh, remesh3D
export settime!

abstract type AbstractField{T, Sf} end

Base.eltype(::AbstractField{T}) where {T} = T
dtype(::AbstractField{T}) where {T} = eltype(T)
sftype(::AbstractField{T, Sf}) where {T, Sf} = Sf
ndofs(field::AbstractField) = length(field.u)
gettags(field::AbstractField) = field.tags

function init!(field::AbstractField, dofmap::AbstractDofMap)
    resize!(field.u, ndofs(dofmap))
    resize!(field.assigned, ndofs(dofmap))
    fill!(field.u, field.u_init)
    fill!(field.assigned, false)

    @assert length(field.u_cache) == length(field.t_cache)
    for i in eachindex(field.u_cache)
        resize!(field.u_cache[i], length(field.u))
        field.u_cache[i] .= field.u
        field.t_cache[i] = field.t[]
    end

    maxid = maximum(shape -> shape.id, eachshape(field))
    resize!(field.dofs, maxid)
    for (shape, dofs) in zip(eachshape(field), dofmap.dofs)
        field.dofs[shape.id] = dofs
    end

    return field
end

function DofMap(field::AbstractField; renumber=true)
    dofmap = DofMap(field.shapefunc, eachshape(field); renumber)
    append!(dofmap.tags, gettags(field))
    return dofmap
end

function settime!(field::AbstractField, t)
    field.t[] = t
    return field
end

settime!(notafield::Any, args...) = notafield

################################
# Domain Field
################################

abstract type AbstractDomainField{T, Sf, M, Pr} <: AbstractField{T, Sf} end

Base.ndims(::AbstractDomainField{T, Sf, M, Pr}) where {T, Sf, M, Pr} = ndims(Pr)
nrefdims(::AbstractDomainField{T, Sf, M, Pr}) where {T, Sf, M, Pr} = nrefdims(Pr)

eachshape(field::AbstractDomainField) = eachelement(field.mesh, field.tags)
eachshape(field::AbstractDomainField, tags) = eachelement(field.mesh, tags)
eachshapeindex(field::AbstractDomainField) = eachelementindex(field.mesh, field.tags)
eachshapeindex(field::AbstractDomainField, tags) = eachelementindex(field.mesh, tags)

function interpolate(field::AbstractDomainField{T, ConstantShapeFunction}, ::Integer) where {T}
    # Elemental fields are constant on an element and therefore not interpolated
    evaluate(field)
end

function evaluate(field::AbstractDomainField{T, ConstantShapeFunction}) where {T}
    proxy = field.proxies[threadid()]
    return field.u[proxy.dofs[1]]
end

function defaultquadrature(field::AbstractDomainField)
    # Pay attention to type stability here if you want type stability
    # for zero allocations in methods using defaultquadrature()
    N = getorder(defaultquadrature(elementtype(field.mesh)))
    M = getorder(defaultquadrature(field.shapefunc))
    return GaussQuadrature(max(N, M))
end

struct Field{T, Sf, M, Pr} <: AbstractDomainField{T, Sf, M, Pr}
    u::Vector{T}
    shapefunc::Sf
    t::RefValue{Float64}
    u_init::T
    u_cache::Vector{Vector{T}}
    t_cache::Vector{Float64}
    tags::Vector{Symbol}
    dofs::Vector{Vector{Int}}
    dofoffset::RefValue{Int}
    mesh::M
    proxies::Vector{Pr}
    assigned::Vector{UInt8}
    outputname::String
end

function Field{T}(mesh::AbstractMesh, shapefunc::AbstractShapeFunctions;
        init=zero(T), cachelength=0, domains=eachdomaintag(mesh),
        outputname="Unnamed Field ($T)") where {T}

    u = T[]
    t = Ref(0.0)
    u_cache = [deepcopy(u) for _ in 1:cachelength]
    t_cache = [0.0 for _ in 1:cachelength]
    dofs = Vector{Int}[]
    dofoffset = Ref(-1)
    proxies = ScalarProxy{ndims(mesh), ndims(mesh), typeof(shapefunc)}[]
    assigned = UInt8[]

    Field{T, typeof(shapefunc), typeof(mesh), eltype(proxies)}(
        u,
        shapefunc,
        t,
        init,
        u_cache,
        t_cache,
        collect(domains),
        dofs,
        dofoffset,
        mesh,
        proxies,
        assigned,
        outputname,
    )
end

const ScalarField{T} = Field{T} where {T <: Number}
const VectorField{N, T} = Field{SVector{N, T}} where {N, T <: Number}

ScalarField(args...; kwargs...) = Field{Float64}(args...; kwargs...)
VectorField{N}(args...; kwargs...) where {N} = Field{SVector{N, Float64}}(args...; kwargs...)

ncomps(::ScalarField) = 1
ncomps(::VectorField{N}) where {N} = N

function Base.show(io::Core.IO, ::MIME"text/plain", f::AbstractDomainField)
    if isempty(f.u)
        println(io, "Uninitialized Field{$(eltype(f)), ...} with:")
    else
        println(io, "Field{$(eltype(f)), ...} with:")
    end
    println(io, "  Interpolation:        $(sftype(f))")
    println(io, "  Number of dofs:       $(ndofs(f))")
    println(io, "  Domains:              $(f.tags)")
    print(io,   "  Output name:          $(f.outputname)")
    if isempty(f.u)
        print(io, "\n  Min, Max:             n/a")
    else
        if eltype(f) <: Real
            print(io, "\n  Min, Max:             $(extrema(f.u))")
        elseif eltype(f) <: Complex
            print(io, "\n  Min, Max (abs):       $(extrema(abs, f.u))")
        elseif eltype(f) <: StaticVector{<:Any, <:Real}
            for i in 1:length(eltype(f))
            print(io, "\n  Min, Max [$i]:        $(extrema(x -> x[i], f.u))")
            end
            print(io, "\n  Min, Max (norm):      $(extrema(norm, f.u))")
        elseif eltype(f) <: StaticVector{<:Any, <:Complex}
            for i in 1:length(eltype(f))
            print(io, "\n  Min, Max [$i] (abs):  $(extrema(x -> abs(x[i]), f.u))")
            end
        else
            print(io, "\n  Min, Max (norm):      not implemented")
        end
    end
end

##############################
# Boundary Field
##############################

abstract type AbstractBoundaryField{T, Sf, M, Pr} <: AbstractField{T, Sf} end

Base.ndims(::AbstractBoundaryField{T, Sf, M, Pr}) where {T, Sf, M, Pr} = ndims(Pr)
nrefdims(::AbstractBoundaryField{T, Sf, M, Pr}) where {T, Sf, M, Pr} = nrefdims(Pr)

eachshape(field::AbstractBoundaryField) = eachfacet(field.mesh, field.tags)
eachshape(field::AbstractBoundaryField, tags) = eachfacet(field.mesh, tags)
eachshapeindex(field::AbstractBoundaryField) = eachfacetindex(field.mesh, field.tags)
eachshapeindex(field::AbstractBoundaryField, tags) = eachfacetindex(field.mesh, tags)

function defaultquadrature(field::AbstractBoundaryField)
    # Pay attention to type stability here if you want type stability
    # for zero allocations in methods using defaultquadrature()
    N = getorder(defaultquadrature(facettype(field.mesh)))
    M = getorder(defaultquadrature(field.shapefunc))
    return GaussQuadrature(max(N, M))
end

# ------------------------------ Scalar Field ------------------------------

struct BoundaryField{T, Sf, M, Pr} <: AbstractBoundaryField{T, Sf, M, Pr}
    u::Vector{T}
    shapefunc::Sf
    t::RefValue{Float64}
    u_init::T
    u_cache::Vector{Vector{T}}
    t_cache::Vector{Float64}
    tags::Vector{Symbol}
    dofs::Vector{Vector{Int}}
    dofoffset::RefValue{Int}
    mesh::M
    proxies::Vector{Pr}
    assigned::Vector{UInt8}
    outputname::String
end

function BoundaryField{T}(mesh::AbstractMesh, shapefunc::AbstractShapeFunctions;
        init=zero(T), cachelength=0, boundaries=eachboundarytag(mesh),
        outputname="Unnamed boundary field ($T)") where {T}

    u = T[]
    t = Ref(0.0)
    u_cache = [deepcopy(u) for _ in 1:cachelength]
    t_cache = [0.0 for _ in 1:cachelength]
    dofs = Vector{Int}[]
    dofoffset = Ref(-1)
    proxies = ScalarProxy{ndims(mesh), ndims(mesh) - 1, typeof(shapefunc)}[]
    assigned = UInt8[]

    BoundaryField{T, typeof(shapefunc), typeof(mesh), eltype(proxies)}(
        u,
        shapefunc,
        t,
        init,
        u_cache,
        t_cache,
        collect(boundaries),
        dofs,
        dofoffset,
        mesh,
        proxies,
        assigned,
        outputname,
    )
end

const ScalarBoundaryField{T} = BoundaryField{T} where {T <: Number}
const VectorBoundaryField{N, T} = BoundaryField{SVector{N, T}} where {N, T <: Number}

ScalarBoundaryField(args...; kwargs...) = BoundaryField{Float64}(args...; kwargs...)
VectorBoundaryField{N}(args...; kwargs...) where {N} = BoundaryField{SVector{N, Float64}}(args...; kwargs...)

ncomps(::ScalarBoundaryField) = 1
ncomps(::VectorBoundaryField{N}) where {N} = N

function Base.show(io::Core.IO, ::MIME"text/plain", f::AbstractBoundaryField)
    if isempty(f.u)
        println(io, "Uninitialized BoundaryField{$(eltype(f)), ...} with:")
    else
        println(io, "BoundaryField{$(eltype(f)), ...} with:")
    end
    println(io, "  Interpolation:        $(sftype(f))")
    println(io, "  Number of dofs:       $(ndofs(f))")
    println(io, "  Boundaries:           $(f.tags)")
    print(io,   "  Output name:          $(f.outputname)")
    if isempty(f.u)
        print(io, "\n  Min, Max:             n/a")
    else
        if eltype(f) <: Real
            print(io, "\n  Min, Max:             $(extrema(f.u))")
        elseif eltype(f) <: Complex
            print(io, "\n  Min, Max (abs):       $(extrema(abs, f.u))")
        elseif eltype(f) <: StaticVector{<:Any, <:Real}
            for i in 1:length(eltype(f))
            print(io, "\n  Min, Max [$i]:        $(extrema(x -> x[i], f.u))")
            end
            print(io, "\n  Min, Max (norm):      $(extrema(norm, f.u))")
        elseif eltype(f) <: StaticVector{<:Any, <:Complex}
            for i in 1:length(eltype(f))
            print(io, "\n  Min, Max [$i] (abs):  $(extrema(x -> abs(x[i]), f.u))")
            end
        else
            print(io, "\n  Min, Max (norm):      not implemented")
        end
    end
end

##############################
# Misc
##############################

"""
"""
eachfield(fields::NamedTuple) = ((name, field) for (name, field) in pairs(fields))

"""
    _convert2isoparametric(::Lagrange{N}, field; renumber=true) -> coordinates, shapes, u

Return `coordinates`, `shapes`, and `u` which represent an isoparametric `Lagrange{N}`
discretization of `field`.

An isoparametric discretization is a discretization where the geometric discretization
and the discretization of `u` are the same (here: Lagrange shape function of the same order).
The dofs of the field are the same as the nodes of the mesh.

Useful for plotting and postprocessing. Most visualizations require that the coordinates
coincide with the dofs of the field.
"""
function _convert2isoparametric(::Lagrange{N}, field::AbstractField; renumber=true) where {N}
    coordinates = field.mesh.coordinates
    shapeindices = eachshapeindex(field)
    shapes = eachshape(field)

    isocoordinates, isoshapes = remesh(Lagrange(N), coordinates, shapes; renumber)
    @assert length(isoshapes) == count(i -> true, shapeindices)

    unassigned = trues(length(isocoordinates))
    u_iso = similar(field.u, length(isocoordinates))
    for (shapeindex, isoshape) in zip(shapeindices, isoshapes)
        for (i, dof) in enumerate(isoshape.n)
            if unassigned[dof]
                ξ = dof_refcoordinates(typeof(isoshape), Lagrange(N), i)
                setproxy_point_evaluation!(field, shapeindex, ξ)
                u_iso[dof] = interpolate(field, 1)
            end
        end
    end

    return isocoordinates, isoshapes, u_iso
end

"""
    _convert2isoparametric1(field) -> newshapes, x_indices, u_indices
    _convert2isoparametric1(shapes, dofs)

Return `newshapes`, `x_indices` and `u_indices` such that
`field.mesh.coordinates[x_indices]`, `newshapes` and `field.u[u_indices]` represent
an isoparametric Lagrange{1} discretization of `field`.

This is useful when a repeated fast conversion of `field` is necessary. In this case,
conversion is only `field.mesh.coordinates[x_indices]` and `field.u[u_indices]`
while `newshapes` are unchanged.

Note: It is assumed that the vertices and vertexdofs always have the same order in the
node and dof ordering scheme for all shapes and Lagrange orders.
"""
function _convert2isoparametric1(field::AbstractField{<:Any, <:Lagrange})
    shapes = eachshape(field)
    dofs = field.dofs
    _convert2isoparametric1(shapes, dofs)
end

function _convert2isoparametric1(shapes, dofs)
    shapes_iso, x_indicies_iso = remesh1(shapes) # Could use remesh1()

    maxindex = maximum(s -> maximum(s.n), shapes_iso)
    u_indices_iso = zeros(Int, maxindex)

    for (shape_iso, shape) in zip(shapes_iso, shapes)
        for i in 1:nvertices(shape_iso)
            u_indices_iso[shape_iso.n[i]] = dofs[shape.id][vertexdofs(typeof(shape), Lagrange(1), i)]
        end
    end

    return shapes_iso, x_indicies_iso, u_indices_iso
end

"""
    _convert2isoparametricfaces1(field) -> newshapes, x_indices, u_indices

TODO
"""
function _convert2isoparametricfaces1(field::AbstractField)
    if nrefdims(field) == 2
        return _convert2isoparametric1(field)
    elseif nrefdims(field) == 3
        mesh = field.mesh
        boundarytags = unique(boundarytag for domaintag in field.tags
            for boundarytag in mesh.domain2boundaries[domaintag])

        facets = [setid(firstorder(facet), i) for (i, facet) in enumerate(eachfacet(mesh, boundarytags))]

        dofs = Vector{Int}[]
        for facetindex in eachfacetindex(mesh, boundarytags)
            localfacetdofs = facetdofs(typeof(getelement(mesh, facetindex)), field.shapefunc, facetindex.ilocal)
            push!(dofs, field.dofs[facetindex.elementid][localfacetdofs])
        end

        return _convert2isoparametric1(facets, dofs)
    else
        error("only 2D and 3D fields are supported.")
    end
end

"""
    _convert2faces1(field) -> newshapes, x_indices, u_indices
"""
_convert2faces1(field::AbstractField) = _convert2isoparametricfaces1(field)

function _convert2faces1(field::AbstractField{<:Any, ConstantShapeFunction})
    dofs = field.dofs

    if nrefdims(field) == 2
        shapes = collect(eachshape(field))
        u_indices = [dofs[shape.id][1] for shape in shapes]
    elseif nrefdims(field) == 3
        mesh = field.mesh
        boundarytags = unique(boundarytag for domaintag in field.tags
            for boundarytag in mesh.domain2boundaries[domaintag])
        shapes = collect(eachfacet(mesh, boundarytags))
        u_indices = [dofs[facetindex.elementid][1] for facetindex in eachfacetindex(mesh, boundarytags)]
    else
        error("only 2D and 3D fields are supported.")
    end

    newshapes, x_indices = remesh1(shapes)

    return newshapes, x_indices, u_indices
end

"""
    remesh(::Lagrange{N}, coordinates, shapes; renumber=false) -> newcoordinates, newshapes

Remesh the mesh defined by `coordinates` and `shapes` using geometric
order `N` (`Lagrange(N)` shape functions) and return the `newcoordinates` and
`newshapes`.

Note: The shape id of the new shapes is rebased to 1, i.e. `newshape.id` corresponds
to its index in vector `newshapes`.

With keyword `renumber=true` the node indices in `newshapes` are renumbered using the
`cuthill_mckee_renumbering` algorithm.

# See also
[`remesh3D`](@ref)

# Example
```julia-repl
julia> coordinates = [SA[0.0], SA[1.0], SA[2.0], SA[3.0]];

julia> shapes = [Edge2(2, 3), Edge2(3, 4)];

julia> newcoordinates, newshapes = remesh(Lagrange(2), coordinates, shapes);

julia> newcoordinates
5-element Vector{SVector{1, Float64}}:
 [0.0]
 [1.0]
 [0.5]
 [2.0]
 [1.5]

julia> newshapes
2-element Vector{Edge3}:
 Edge3([1, 2, 3], 1)
 Edge3([2, 4, 5], 2)
```
"""
function remesh(::Lagrange{N}, coordinates, shapes; renumber=false) where N
    dofmap = DofMap(Lagrange(N), shapes; renumber)
    newshapes = [shapetype(typeof(shape), Lagrange(N))(dofmap.dofs[i]; id=i) for (i, shape) in enumerate(shapes)]

    unassigned = trues(ndofs(dofmap))
    newcoordinates = Vector{eltype(coordinates)}(undef, ndofs(dofmap))

    for (shape, newshape) in zip(shapes, newshapes)
        for (i, node) in enumerate(newshape.n)
            if unassigned[node]
                ξ = dof_refcoordinates(typeof(newshape), Lagrange(N), i)
                shapecoordinates = @view coordinates[shape.n]
                x = to_globalcoordinates(ξ, typeof(shape), shapecoordinates)
                newcoordinates[node] = x
                unassigned[node] = false
            end
        end
    end

    return newcoordinates, newshapes
end

"""
    remesh3D(::Lagrange{N}, coordinates, shapes; renumber=false) -> newcoordinates3D, newshapes

Similar to [`remesh`](@ref), but the newcoordinates are always padded with zeros to 3D coordinates.

# Example
```julia-repl
julia> coordinates = [SA[0.0], SA[1.0], SA[2.0], SA[3.0]];

julia> shapes = [Edge2(2, 3), Edge2(3, 4)];

julia> newcoordinates, newshapes = remesh3D(Lagrange(2), coordinates, shapes);

julia> newcoordinates
5-element Vector{SVector{3, Float64}}:
 [1.0, 0.0, 0.0]
 [2.0, 0.0, 0.0]
 [1.5, 0.0, 0.0]
 [3.0, 0.0, 0.0]
 [2.5, 0.0, 0.0]

julia> newshapes
2-element Vector{Edge3}:
 Edge3([1, 2, 3], 1)
 Edge3([2, 4, 5], 2)
```
"""
function remesh3D(::Lagrange{N}, coordinates, shapes; renumber=false) where N
    dofmap = DofMap(Lagrange(N), shapes; renumber)
    newshapes = [shapetype(typeof(shape), Lagrange(N))(dofmap.dofs[i]; id=i) for (i, shape) in enumerate(shapes)]

    T = typeof(pad3D(first(coordinates)))
    unassigned = trues(ndofs(dofmap))
    newcoordinates = Vector{T}(undef, ndofs(dofmap))

    for (shape, newshape) in zip(shapes, newshapes)
        for (i, node) in enumerate(newshape.n)
            if unassigned[node]
                ξ = dof_refcoordinates(typeof(newshape), Lagrange(N), i)
                shapecoordinates = @view coordinates[shape.n]
                x = to_globalcoordinates(ξ, typeof(shape), shapecoordinates)
                newcoordinates[node] = pad3D(x)
                unassigned[node] = false
            end
        end
    end

    return newcoordinates, newshapes
end

"""
    remesh1(shapes; renumber=false) -> newshapes, indices

Return `newshapes` and `indices` such that `newshapes` and `coordinates[indices]`
are a geometrically first order mesh (`Lagrange{1}`) based on the mesh defined by `shapes`
`coordinates`.

This is useful when a repeated fast mesh conversion is necessary. In this case,
conversion is only `coordinates[indices]` while `newshapes` are unchanged.

Note: It is assumed that the vertices and vertexdofs always have the same order in the
node and dof ordering scheme for all shapes and Lagrange orders.

Note: The shape id of the new shapes is rebased to 1, i.e. `newshape.id` corresponds
to its index in vector `newshapes`.

With keyword `renumber=true` the node indices in `newshapes` are renumbered using the
`cuthill_mckee_renumbering` algorithm.

# See also
[`remesh`](@ref)

# Example
```julia-repl
julia> coordinates = [SA[0.0], SA[1.0], SA[0.5], SA[2.0], SA[1.5]];

julia> shapes = [Edge3([1, 2, 3]), Edge3([2, 4, 5])];

julia> newshapes, indices = remesh1(shapes);

julia> newshapes
2-element Vector{Edge2}:
 Edge2([1, 2], 1)
 Edge2([2, 3], 2)

julia> (coordinates[indices], newshapes) == remesh(Lagrange(1), coordinates, shapes)
 true
```
"""
function remesh1(shapes; renumber=false)
    dofmap = DofMap(Lagrange(1), shapes; renumber)
    newshapes = [shapetype(typeof(shape), Lagrange(1))(dofmap.dofs[i]; id=i) for (i, shape) in enumerate(shapes)]

    unassigned = trues(ndofs(dofmap))
    indices = zeros(Int, ndofs(dofmap))

    for (shape, newshape) in zip(shapes, newshapes)
        for (i, node) in enumerate(newshape.n)
            if unassigned[node]
                indices[node] = shape.n[i]
                unassigned[node] = false
            end
        end
    end

    return newshapes, indices
end

##############################
# Field Handler
##############################

abstract type AbstractFieldHandler end

"""
"""
struct FieldHandler{F} <: AbstractFieldHandler
    fields::F
end

function FieldHandler!(fields; renumber=true, verbose=false)
    for (name, field) in pairs(fields)
        field isa AbstractField || continue

        for (name2, field2) in pairs(fields)
            field2 isa AbstractField || continue
            if field !== field2 && !isempty(field2.proxies) && is_similar_proxy(field, field2)
                verbose && println("Reusing proxy of field `$name2` for field `$name`.")
                resize!(field.proxies, length(field2.proxies))
                field.proxies .= field2.proxies
                break
            end
            append!(field.proxies, eltype(field.proxies)(field.shapefunc) for _ in 1:nthreads())
        end

        for (name2, field2) in pairs(fields)
            field2 isa AbstractField || continue
            if field !== field2 && !isempty(field2.dofs) && is_similar_dofmap(field, field2)
                verbose && println("Reusing dofmap of field `$name2` for field `$name`.")
                resize!(field.dofs, length(field2.dofs))
                field.dofs .= field2.dofs
                resize!(field.u, length(field2.u))
                fill!(field.u, field.u_init)
                resize!(field.assigned, length(field2.assigned))
                fill!(field.assigned, false)
                break
            end
            new_dofmap = DofMap(field; renumber)
            init!(field, new_dofmap)
        end
    end

    return FieldHandler{typeof(fields)}(fields)
end

FieldHandler!(field::AbstractField, args...; kwargs...) = FieldHandler!((u = field,), args...; kwargs...)

function Base.show(io::Core.IO, ::MIME"text/plain", fh::FieldHandler)
    println(io, "FieldHandler handling $(length(fh.fields)) fields with:")
    println(io, "  Field names:           $(propertynames(fh.fields))")
end

function is_similar_proxy(field1, field2)
    if eltype(field1.proxies) === eltype(field2.proxies)
        return true
    end
    return false
end

function is_similar_dofmap(field1, field2)
    if field1.mesh === field2.mesh
        if field1.shapefunc === field2.shapefunc
            if field1.tags == field2.tags
                return true
            end
        end
    end
    return false
end

function defaultquadrature(fh::AbstractFieldHandler)
    fields = (field for field in fh.fields if field isa AbstractField)
    N = maximum(getorder(defaultquadrature(field)) for field in fields)
    return GaussQuadrature(N)
end

settime!(fh::AbstractFieldHandler, args...) = settime!(fh.fields, args...)

#####################################
# Set Proxy Integration
#####################################

function setproxy_integration!(field::AbstractDomainField, index::ElementIndex, quadrature=defaultquadrature(field))
    proxy = field.proxies[threadid()]
    if proxy._mode[] == 0 && proxy._quadorder[] == getorder(quadrature) && proxy._shapeid[] == index.shapeid
        return field
    end
    shape = field.mesh[index]
    _init_element_integration!(proxy, typeof(shape), quadrature)
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    _set_element_integration!(proxy, typeof(shape), shapecoordinates, quadrature)
    setdofs!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

function setproxy_integration!(field::AbstractDomainField, index::FacetIndex, quadrature=defaultquadrature(field))
    proxy = field.proxies[threadid()]
    if proxy._mode[] == index.ilocal && proxy._quadorder[] == getorder(quadrature) && proxy._shapeid[] == index.elementid
        return field
    end
    shape = getelement(field.mesh, index)
    _init_facet_integration!(proxy, typeof(shape), quadrature)
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    _set_facet_integration!(proxy, typeof(shape), shapecoordinates, index.ilocal, quadrature)
    setdofs!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

function setproxy_integration!(field::AbstractBoundaryField, index::FacetIndex, quadrature=defaultquadrature(field))
    proxy = field.proxies[threadid()]
    if proxy._mode[] == 0 && proxy._quadorder[] == getorder(quadrature) && proxy._shapeid[] == index.shapeid
        return field
    end
    shape = field.mesh[index]
    _init_element_integration!(proxy, typeof(shape), quadrature)
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    _set_element_integration!(proxy, typeof(shape), shapecoordinates, quadrature)
    setdofs!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

function setproxy_integration!(notfield::Any, args...)
    return notfield
end

#####################################
# Set Proxy Dof Evaluation
#####################################

function setproxy_dof_evaluation!(field::AbstractDomainField, index::ElementIndex)
    proxy = field.proxies[threadid()]
    if proxy._mode[] == 10 && proxy._shapeid[] == index.shapeid
        return field
    end
    shape = field.mesh[index]
    _init_element_dof_evaluation!(proxy, typeof(shape))
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    _set_element_dof_evaluation!(proxy, typeof(shape), shapecoordinates)
    setdofs_noperm!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

function setproxy_dof_evaluation!(field::AbstractDomainField, index::FacetIndex)
    proxy = field.proxies[threadid()]
    if proxy._mode[] == -index.ilocal && proxy._shapeid[] == index.elementid
        return field
    end
    shape = getelement(field.mesh, index)
    _init_facet_dof_evaluation!(proxy, typeof(shape))
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    _set_facet_dof_evaluation!(proxy, typeof(shape), shapecoordinates, index.ilocal)
    setdofs_noperm!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

function setproxy_dof_evaluation!(field::AbstractBoundaryField, index::FacetIndex)
    proxy = field.proxies[threadid()]
    if proxy._mode[] == 10 && proxy._shapeid[] == index.shapeid
        return field
    end
    shape = field.mesh[index]
    _init_element_dof_evaluation!(proxy, typeof(shape))
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    _set_element_dof_evaluation!(proxy, typeof(shape), shapecoordinates)
    setdofs_noperm!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

function setproxy_dof_evaluation!(notfield::Any, args...)
    return notfield
end

#####################################
# Set Proxy Multi Point Evaluation
#####################################

function setproxy_point_evaluation!(notfield::Any, args...)
    return notfield
end

function setproxy_point_evaluation!(field::AbstractDomainField, index::ElementIndex)
    proxy = field.proxies[threadid()]
    if proxy._mode[] == 20 && proxy._shapeid[] == index.shapeid
        return field
    end
    shape = field.mesh[index]
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    set_point_evaluation!(proxy, typeof(shape), shapecoordinates)
    setdofs_noperm!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

function setproxy_point_evaluation!(field::AbstractBoundaryField, index::FacetIndex)
    proxy = field.proxies[threadid()]
    if proxy._mode[] == 20 && proxy._shapeid[] == index.shapeid
        return field
    end
    shape = field.mesh[index]
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    set_point_evaluation!(proxy, typeof(shape), shapecoordinates)
    setdofs_noperm!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

function initproxy_point_evaluation!(field::AbstractField, index::AbstractMeshIndex, ξs::AbstractVector{<:AbstractVector})
    shape = field.mesh[index]
    initproxy_point_evaluation!(field, typeof(shape), ξs)
end

function initproxy_point_evaluation!(field::AbstractField, ::Type{S}, ξs::AbstractVector{<:AbstractVector}) where {S<:AbstractShape}
    proxy = field.proxies[threadid()]
    init_point_evaluation!(proxy, S, ξs)
    return field
end

#####################################
# Set Proxy Single Point Evaluation
#####################################

function setproxy_point_evaluation!(field::AbstractDomainField, index::ElementIndex, ξ::AbstractVector{<:Number})
    proxy = field.proxies[threadid()]
    if proxy._mode[] == 21 && proxy._shapeid[] == index.shapeid && proxy._evalpoints[1][1] == ξ
        return field
    end
    shape = field.mesh[index]
    _init_point_evaluation!(proxy, typeof(shape))
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    _set_point_evaluation!(proxy, typeof(shape), shapecoordinates, ξ)
    setdofs_noperm!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

function setproxy_point_evaluation!(field::AbstractBoundaryField, index::FacetIndex, ξ::AbstractVector{<:Number})
    proxy = field.proxies[threadid()]
    if proxy._mode[] == 21 && proxy._shapeid[] == index.shapeid && proxy._evalpoints[1][1] == ξ
        return field
    end
    shape = field.mesh[index]
    _init_point_evaluation!(proxy, typeof(shape))
    shapecoordinates = @view field.mesh.coordinates[shape.n]
    _set_point_evaluation!(proxy, typeof(shape), shapecoordinates, ξ)
    setdofs_noperm!(proxy, field.dofs[shape.id])
    proxy._shapeid[] = shape.id
    return field
end

#####################################
# setproxy Dispatch for field tuples
#####################################
# TODO Generated functions

function setproxy_integration!(fields::TupleOfLength{1, Any}, args...)
    f1, = fields
    setproxy_integration!(f1, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{2, Any}, args...)
    f1, f2 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{3, Any}, args...)
    f1, f2, f3 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{4, Any}, args...)
    f1, f2, f3, f4 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{5, Any}, args...)
    f1, f2, f3, f4, f5 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{6, Any}, args...)
    f1, f2, f3, f4, f5, f6 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{7, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{8, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    setproxy_integration!(f8, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{9, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    setproxy_integration!(f8, args...)
    setproxy_integration!(f9, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{10, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    setproxy_integration!(f8, args...)
    setproxy_integration!(f9, args...)
    setproxy_integration!(f10, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{11, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    setproxy_integration!(f8, args...)
    setproxy_integration!(f9, args...)
    setproxy_integration!(f10, args...)
    setproxy_integration!(f11, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{12, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    setproxy_integration!(f8, args...)
    setproxy_integration!(f9, args...)
    setproxy_integration!(f10, args...)
    setproxy_integration!(f11, args...)
    setproxy_integration!(f12, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{13, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    setproxy_integration!(f8, args...)
    setproxy_integration!(f9, args...)
    setproxy_integration!(f10, args...)
    setproxy_integration!(f11, args...)
    setproxy_integration!(f12, args...)
    setproxy_integration!(f13, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{14, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    setproxy_integration!(f8, args...)
    setproxy_integration!(f9, args...)
    setproxy_integration!(f10, args...)
    setproxy_integration!(f11, args...)
    setproxy_integration!(f12, args...)
    setproxy_integration!(f13, args...)
    setproxy_integration!(f14, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{15, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    setproxy_integration!(f8, args...)
    setproxy_integration!(f9, args...)
    setproxy_integration!(f10, args...)
    setproxy_integration!(f11, args...)
    setproxy_integration!(f12, args...)
    setproxy_integration!(f13, args...)
    setproxy_integration!(f14, args...)
    setproxy_integration!(f15, args...)
    return fields
end

function setproxy_integration!(fields::TupleOfLength{16, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16 = fields
    setproxy_integration!(f1, args...)
    setproxy_integration!(f2, args...)
    setproxy_integration!(f3, args...)
    setproxy_integration!(f4, args...)
    setproxy_integration!(f5, args...)
    setproxy_integration!(f6, args...)
    setproxy_integration!(f7, args...)
    setproxy_integration!(f8, args...)
    setproxy_integration!(f9, args...)
    setproxy_integration!(f10, args...)
    setproxy_integration!(f11, args...)
    setproxy_integration!(f12, args...)
    setproxy_integration!(f13, args...)
    setproxy_integration!(f14, args...)
    setproxy_integration!(f15, args...)
    setproxy_integration!(f16, args...)
    return fields
end

# -----------------------------

function setproxy_dof_evaluation!(fields::TupleOfLength{1, Any}, args...)
    f1, = fields
    setproxy_dof_evaluation!(f1, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{2, Any}, args...)
    f1, f2 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{3, Any}, args...)
    f1, f2, f3 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{4, Any}, args...)
    f1, f2, f3, f4 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{5, Any}, args...)
    f1, f2, f3, f4, f5 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{6, Any}, args...)
    f1, f2, f3, f4, f5, f6 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{7, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{8, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    setproxy_dof_evaluation!(f8, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{9, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    setproxy_dof_evaluation!(f8, args...)
    setproxy_dof_evaluation!(f9, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{10, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    setproxy_dof_evaluation!(f8, args...)
    setproxy_dof_evaluation!(f9, args...)
    setproxy_dof_evaluation!(f10, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{11, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    setproxy_dof_evaluation!(f8, args...)
    setproxy_dof_evaluation!(f9, args...)
    setproxy_dof_evaluation!(f10, args...)
    setproxy_dof_evaluation!(f11, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{12, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    setproxy_dof_evaluation!(f8, args...)
    setproxy_dof_evaluation!(f9, args...)
    setproxy_dof_evaluation!(f10, args...)
    setproxy_dof_evaluation!(f11, args...)
    setproxy_dof_evaluation!(f12, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{13, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    setproxy_dof_evaluation!(f8, args...)
    setproxy_dof_evaluation!(f9, args...)
    setproxy_dof_evaluation!(f10, args...)
    setproxy_dof_evaluation!(f11, args...)
    setproxy_dof_evaluation!(f12, args...)
    setproxy_dof_evaluation!(f13, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{14, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    setproxy_dof_evaluation!(f8, args...)
    setproxy_dof_evaluation!(f9, args...)
    setproxy_dof_evaluation!(f10, args...)
    setproxy_dof_evaluation!(f11, args...)
    setproxy_dof_evaluation!(f12, args...)
    setproxy_dof_evaluation!(f13, args...)
    setproxy_dof_evaluation!(f14, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{15, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    setproxy_dof_evaluation!(f8, args...)
    setproxy_dof_evaluation!(f9, args...)
    setproxy_dof_evaluation!(f10, args...)
    setproxy_dof_evaluation!(f11, args...)
    setproxy_dof_evaluation!(f12, args...)
    setproxy_dof_evaluation!(f13, args...)
    setproxy_dof_evaluation!(f14, args...)
    setproxy_dof_evaluation!(f15, args...)
    return fields
end

function setproxy_dof_evaluation!(fields::TupleOfLength{16, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16 = fields
    setproxy_dof_evaluation!(f1, args...)
    setproxy_dof_evaluation!(f2, args...)
    setproxy_dof_evaluation!(f3, args...)
    setproxy_dof_evaluation!(f4, args...)
    setproxy_dof_evaluation!(f5, args...)
    setproxy_dof_evaluation!(f6, args...)
    setproxy_dof_evaluation!(f7, args...)
    setproxy_dof_evaluation!(f8, args...)
    setproxy_dof_evaluation!(f9, args...)
    setproxy_dof_evaluation!(f10, args...)
    setproxy_dof_evaluation!(f11, args...)
    setproxy_dof_evaluation!(f12, args...)
    setproxy_dof_evaluation!(f13, args...)
    setproxy_dof_evaluation!(f14, args...)
    setproxy_dof_evaluation!(f15, args...)
    setproxy_dof_evaluation!(f16, args...)
    return fields
end

# -----------------------------

function setproxy_point_evaluation!(fields::TupleOfLength{1, Any}, args...)
    f1, = fields
    setproxy_point_evaluation!(f1, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{2, Any}, args...)
    f1, f2 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{3, Any}, args...)
    f1, f2, f3 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{4, Any}, args...)
    f1, f2, f3, f4 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{5, Any}, args...)
    f1, f2, f3, f4, f5 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{6, Any}, args...)
    f1, f2, f3, f4, f5, f6 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{7, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{8, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    setproxy_point_evaluation!(f8, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{9, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    setproxy_point_evaluation!(f8, args...)
    setproxy_point_evaluation!(f9, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{10, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    setproxy_point_evaluation!(f8, args...)
    setproxy_point_evaluation!(f9, args...)
    setproxy_point_evaluation!(f10, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{11, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    setproxy_point_evaluation!(f8, args...)
    setproxy_point_evaluation!(f9, args...)
    setproxy_point_evaluation!(f10, args...)
    setproxy_point_evaluation!(f11, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{12, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    setproxy_point_evaluation!(f8, args...)
    setproxy_point_evaluation!(f9, args...)
    setproxy_point_evaluation!(f10, args...)
    setproxy_point_evaluation!(f11, args...)
    setproxy_point_evaluation!(f12, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{13, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    setproxy_point_evaluation!(f8, args...)
    setproxy_point_evaluation!(f9, args...)
    setproxy_point_evaluation!(f10, args...)
    setproxy_point_evaluation!(f11, args...)
    setproxy_point_evaluation!(f12, args...)
    setproxy_point_evaluation!(f13, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{14, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    setproxy_point_evaluation!(f8, args...)
    setproxy_point_evaluation!(f9, args...)
    setproxy_point_evaluation!(f10, args...)
    setproxy_point_evaluation!(f11, args...)
    setproxy_point_evaluation!(f12, args...)
    setproxy_point_evaluation!(f13, args...)
    setproxy_point_evaluation!(f14, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{15, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    setproxy_point_evaluation!(f8, args...)
    setproxy_point_evaluation!(f9, args...)
    setproxy_point_evaluation!(f10, args...)
    setproxy_point_evaluation!(f11, args...)
    setproxy_point_evaluation!(f12, args...)
    setproxy_point_evaluation!(f13, args...)
    setproxy_point_evaluation!(f14, args...)
    setproxy_point_evaluation!(f15, args...)
    return fields
end

function setproxy_point_evaluation!(fields::TupleOfLength{16, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16 = fields
    setproxy_point_evaluation!(f1, args...)
    setproxy_point_evaluation!(f2, args...)
    setproxy_point_evaluation!(f3, args...)
    setproxy_point_evaluation!(f4, args...)
    setproxy_point_evaluation!(f5, args...)
    setproxy_point_evaluation!(f6, args...)
    setproxy_point_evaluation!(f7, args...)
    setproxy_point_evaluation!(f8, args...)
    setproxy_point_evaluation!(f9, args...)
    setproxy_point_evaluation!(f10, args...)
    setproxy_point_evaluation!(f11, args...)
    setproxy_point_evaluation!(f12, args...)
    setproxy_point_evaluation!(f13, args...)
    setproxy_point_evaluation!(f14, args...)
    setproxy_point_evaluation!(f15, args...)
    setproxy_point_evaluation!(f16, args...)
    return fields
end

# -----------------------------

function initproxy_point_evaluation!(fields::TupleOfLength{1, Any}, args...)
    f1, = fields
    initproxy_point_evaluation!(f1, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{2, Any}, args...)
    f1, f2 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{3, Any}, args...)
    f1, f2, f3 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{4, Any}, args...)
    f1, f2, f3, f4 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{5, Any}, args...)
    f1, f2, f3, f4, f5 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{6, Any}, args...)
    f1, f2, f3, f4, f5, f6 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{7, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{8, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    initproxy_point_evaluation!(f8, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{9, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    initproxy_point_evaluation!(f8, args...)
    initproxy_point_evaluation!(f9, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{10, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    initproxy_point_evaluation!(f8, args...)
    initproxy_point_evaluation!(f9, args...)
    initproxy_point_evaluation!(f10, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{11, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    initproxy_point_evaluation!(f8, args...)
    initproxy_point_evaluation!(f9, args...)
    initproxy_point_evaluation!(f10, args...)
    initproxy_point_evaluation!(f11, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{12, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    initproxy_point_evaluation!(f8, args...)
    initproxy_point_evaluation!(f9, args...)
    initproxy_point_evaluation!(f10, args...)
    initproxy_point_evaluation!(f11, args...)
    initproxy_point_evaluation!(f12, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{13, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    initproxy_point_evaluation!(f8, args...)
    initproxy_point_evaluation!(f9, args...)
    initproxy_point_evaluation!(f10, args...)
    initproxy_point_evaluation!(f11, args...)
    initproxy_point_evaluation!(f12, args...)
    initproxy_point_evaluation!(f13, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{14, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    initproxy_point_evaluation!(f8, args...)
    initproxy_point_evaluation!(f9, args...)
    initproxy_point_evaluation!(f10, args...)
    initproxy_point_evaluation!(f11, args...)
    initproxy_point_evaluation!(f12, args...)
    initproxy_point_evaluation!(f13, args...)
    initproxy_point_evaluation!(f14, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{15, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    initproxy_point_evaluation!(f8, args...)
    initproxy_point_evaluation!(f9, args...)
    initproxy_point_evaluation!(f10, args...)
    initproxy_point_evaluation!(f11, args...)
    initproxy_point_evaluation!(f12, args...)
    initproxy_point_evaluation!(f13, args...)
    initproxy_point_evaluation!(f14, args...)
    initproxy_point_evaluation!(f15, args...)
    return fields
end

function initproxy_point_evaluation!(fields::TupleOfLength{16, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16 = fields
    initproxy_point_evaluation!(f1, args...)
    initproxy_point_evaluation!(f2, args...)
    initproxy_point_evaluation!(f3, args...)
    initproxy_point_evaluation!(f4, args...)
    initproxy_point_evaluation!(f5, args...)
    initproxy_point_evaluation!(f6, args...)
    initproxy_point_evaluation!(f7, args...)
    initproxy_point_evaluation!(f8, args...)
    initproxy_point_evaluation!(f9, args...)
    initproxy_point_evaluation!(f10, args...)
    initproxy_point_evaluation!(f11, args...)
    initproxy_point_evaluation!(f12, args...)
    initproxy_point_evaluation!(f13, args...)
    initproxy_point_evaluation!(f14, args...)
    initproxy_point_evaluation!(f15, args...)
    initproxy_point_evaluation!(f16, args...)
    return fields
end

# -----------------------------

function setproxy!(fields::TupleOfLength{1, Any}, args...)
    f1, = fields
    setproxy!(f1, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{2, Any}, args...)
    f1, f2 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{3, Any}, args...)
    f1, f2, f3 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{4, Any}, args...)
    f1, f2, f3, f4 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{5, Any}, args...)
    f1, f2, f3, f4, f5 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{6, Any}, args...)
    f1, f2, f3, f4, f5, f6 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{7, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{8, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    setproxy!(f8, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{9, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    setproxy!(f8, args...)
    setproxy!(f9, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{10, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    setproxy!(f8, args...)
    setproxy!(f9, args...)
    setproxy!(f10, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{11, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    setproxy!(f8, args...)
    setproxy!(f9, args...)
    setproxy!(f10, args...)
    setproxy!(f11, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{12, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    setproxy!(f8, args...)
    setproxy!(f9, args...)
    setproxy!(f10, args...)
    setproxy!(f11, args...)
    setproxy!(f12, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{13, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    setproxy!(f8, args...)
    setproxy!(f9, args...)
    setproxy!(f10, args...)
    setproxy!(f11, args...)
    setproxy!(f12, args...)
    setproxy!(f13, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{14, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    setproxy!(f8, args...)
    setproxy!(f9, args...)
    setproxy!(f10, args...)
    setproxy!(f11, args...)
    setproxy!(f12, args...)
    setproxy!(f13, args...)
    setproxy!(f14, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{15, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    setproxy!(f8, args...)
    setproxy!(f9, args...)
    setproxy!(f10, args...)
    setproxy!(f11, args...)
    setproxy!(f12, args...)
    setproxy!(f13, args...)
    setproxy!(f14, args...)
    setproxy!(f15, args...)
    return fields
end

function setproxy!(fields::TupleOfLength{16, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16 = fields
    setproxy!(f1, args...)
    setproxy!(f2, args...)
    setproxy!(f3, args...)
    setproxy!(f4, args...)
    setproxy!(f5, args...)
    setproxy!(f6, args...)
    setproxy!(f7, args...)
    setproxy!(f8, args...)
    setproxy!(f9, args...)
    setproxy!(f10, args...)
    setproxy!(f11, args...)
    setproxy!(f12, args...)
    setproxy!(f13, args...)
    setproxy!(f14, args...)
    setproxy!(f15, args...)
    setproxy!(f16, args...)
    return fields
end

#####################################
# Dispatch for settime!
#####################################

function settime!(fields::TupleOfLength{1, Any}, args...)
    f1, = fields
    settime!(f1, args...)
    return fields
end

function settime!(fields::TupleOfLength{2, Any}, args...)
    f1, f2 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    return fields
end

function settime!(fields::TupleOfLength{3, Any}, args...)
    f1, f2, f3 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    return fields
end

function settime!(fields::TupleOfLength{4, Any}, args...)
    f1, f2, f3, f4 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    return fields
end

function settime!(fields::TupleOfLength{5, Any}, args...)
    f1, f2, f3, f4, f5 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    return fields
end

function settime!(fields::TupleOfLength{6, Any}, args...)
    f1, f2, f3, f4, f5, f6 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    return fields
end

function settime!(fields::TupleOfLength{7, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    return fields
end

function settime!(fields::TupleOfLength{8, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    settime!(f8, args...)
    return fields
end

function settime!(fields::TupleOfLength{9, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    settime!(f8, args...)
    settime!(f9, args...)
    return fields
end

function settime!(fields::TupleOfLength{10, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    settime!(f8, args...)
    settime!(f9, args...)
    settime!(f10, args...)
    return fields
end

function settime!(fields::TupleOfLength{11, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    settime!(f8, args...)
    settime!(f9, args...)
    settime!(f10, args...)
    settime!(f11, args...)
    return fields
end

function settime!(fields::TupleOfLength{12, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    settime!(f8, args...)
    settime!(f9, args...)
    settime!(f10, args...)
    settime!(f11, args...)
    settime!(f12, args...)
    return fields
end

function settime!(fields::TupleOfLength{13, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    settime!(f8, args...)
    settime!(f9, args...)
    settime!(f10, args...)
    settime!(f11, args...)
    settime!(f12, args...)
    settime!(f13, args...)
    return fields
end

function settime!(fields::TupleOfLength{14, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    settime!(f8, args...)
    settime!(f9, args...)
    settime!(f10, args...)
    settime!(f11, args...)
    settime!(f12, args...)
    settime!(f13, args...)
    settime!(f14, args...)
    return fields
end

function settime!(fields::TupleOfLength{15, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    settime!(f8, args...)
    settime!(f9, args...)
    settime!(f10, args...)
    settime!(f11, args...)
    settime!(f12, args...)
    settime!(f13, args...)
    settime!(f14, args...)
    settime!(f15, args...)
    return fields
end

function settime!(fields::TupleOfLength{16, Any}, args...)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16 = fields
    settime!(f1, args...)
    settime!(f2, args...)
    settime!(f3, args...)
    settime!(f4, args...)
    settime!(f5, args...)
    settime!(f6, args...)
    settime!(f7, args...)
    settime!(f8, args...)
    settime!(f9, args...)
    settime!(f10, args...)
    settime!(f11, args...)
    settime!(f12, args...)
    settime!(f13, args...)
    settime!(f14, args...)
    settime!(f15, args...)
    settime!(f16, args...)
    return fields
end
