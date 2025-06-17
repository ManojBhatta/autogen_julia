export VTKWriter
export vtkcelltype
export enable!
export disable!
export reset_counter!
export write!
export write_last!
export write_forced!
export create_vtk_grid
export writevtk
export outputvars
export outputnames
export outputvars_element
export outputnames_element

#############################
# Simple writevtk
#############################

"""
    writevtk(filename, field::AbstractField)

Write `field` to a VTK file `filename`.

When writing multiple fields, solutions of PDEs or when writing
in a loop consider using `VTKWriter` instead which can cache and reuse dofs
and probes for better performance and allows for writing VTK files
as a timeseries.
"""
function writevtk(filename, field::AbstractField; order=defaultvtkorder(field))
    writevtk(filename, field, Lagrange(order))
end

function writevtk(filename, field::AbstractField, shapefunc::Lagrange)
    vtkpoints, vtkcells, u = convert2vtkvars(shapefunc, field)
    mkpath(dirname(filename))
    vtk = vtk_grid(filename, vtkpoints, vtkcells)
    vtk[field.outputname] = u
    vtk_save(vtk)
    return vtk
end

function writevtk(filename, mesh::AbstractMesh; domains=eachdomaintag(mesh))
    writevtk(filename, mesh.coordinates, eachelement(mesh, domains))
end


"""
    writevtk(filename, coordinates, shapes)
    writevtk(filename, coordinates, shapes, "fieldname1" => u1, ....)

Discretization must be isoparametric, i.e. the nodes of a shape index both the coordinates
and the dofs of the field `u1`), or `u1` must be cell field, i.e. `u` must be a vector of the
same length as `shapes`. `u1` can also be a `<:Number`.
"""
function writevtk(filename, coordinates, shapes, args...)
    vtk = _writevtk(filename, coordinates, shapes, args)
    vtk_save(vtk)
end

function _writevtk(filename, coordinates, shapes, fieldpairs::Tuple{Vararg{Pair{String, <:Any}}})
    vtkcells = [MeshCell(vtkcelltype(typeof(shape)), shape.n) for shape in shapes]
    mkpath(dirname(filename))
    vtk = vtk_grid(filename, map(pad3D, coordinates), vtkcells)
    for fieldpair in fieldpairs
        vtk[fieldpair.first] = fieldpair.second
    end
    return vtk
end

function convert2vtkvars(sf::Lagrange, field::AbstractField)
    coordinates, shapes, u = _convert2isoparametric(sf, field, renumber=false)
    vtkcells = [MeshCell(vtkcelltype(typeof(shape)), shape.n) for shape in shapes]
    vtkpoints = map(pad3D, coordinates)
    return vtkpoints, vtkcells, u
end

function convert2vtkvars(sf::Lagrange, field::AbstractField{<:Any, ConstantShapeFunction})
    coordinates = field.mesh.coordinates
    shapes = eachshape(field)
    vtkpoints, newshapes = remesh3D(sf, coordinates, shapes, renumber=false)
    vtkcells = [MeshCell(vtkcelltype(typeof(shape)), shape.n) for shape in newshapes]
    return vtkpoints, vtkcells, field.u
end

defaultvtkorder(field::AbstractField) = defaultvtkorder(eachshape(field), field.shapefunc)

function defaultvtkorder(shapes, shapefunc::AbstractShapeFunctions)
    n_geo = getorder(sftype(typeof(first(shapes))))
    if n_geo < 2 && getorder(shapefunc) < 2
        return 1
    else
        return 2
    end
end

#############################
# AbstractWriter
#############################

abstract type AbstractWriter end

function enable!(w::AbstractWriter)
    w.disabled[] = false
end

function disable!(w::AbstractWriter)
    w.disabled[] = true
end

Base.length(::AbstractWriter) = 1
Base.iterate(w::AbstractWriter) = (w, nothing)
Base.iterate(::AbstractWriter, ::Any) = nothing

nframes(w) = length(w.timestamps)

# write!(writer::AbstractWriter, args...; kwargs...) = write!((writer,), args...; kwargs...)

function write!(writers::Tuple{Vararg{AbstractWriter, 1}}, args...; kwargs...)
    writer1, = writers
    write!(writer1, args...; kwargs...)
    return
end

function write!(writers::Tuple{Vararg{AbstractWriter, 2}}, args...; kwargs...)
    writer1, writer2 = writers
    write!(writer1, args...; kwargs...)
    write!(writer2, args...; kwargs...)
    return
end

function write!(writers::Tuple{Vararg{AbstractWriter, 3}}, args...; kwargs...)
    writer1, writer2, writer3 = writers
    write!(writer1, args...; kwargs...)
    write!(writer2, args...; kwargs...)
    write!(writer3, args...; kwargs...)
    return
end

function write!(writers::Tuple{Vararg{AbstractWriter, 4}}, args...; kwargs...)
    writer1, writer2, writer3, writer4 = writers
    write!(writer1, args...; kwargs...)
    write!(writer2, args...; kwargs...)
    write!(writer3, args...; kwargs...)
    write!(writer4, args...; kwargs...)
    return
end

function write!(writers::Tuple{Vararg{AbstractWriter, 5}}, args...; kwargs...)
    writer1, writer2, writer3, writer4, writer5 = writers
    write!(writer1, args...; kwargs...)
    write!(writer2, args...; kwargs...)
    write!(writer3, args...; kwargs...)
    write!(writer4, args...; kwargs...)
    write!(writer5, args...; kwargs...)
    return
end

function write!(writers::Tuple{Vararg{AbstractWriter, 6}}, args...; kwargs...)
    writer1, writer2, writer3, writer4, writer5, writer6 = writers
    write!(writer1, args...; kwargs...)
    write!(writer2, args...; kwargs...)
    write!(writer3, args...; kwargs...)
    write!(writer4, args...; kwargs...)
    write!(writer5, args...; kwargs...)
    write!(writer6, args...; kwargs...)
    return
end

write!(::Nothing, args...) = nothing

"""
"""
reset_counter!(writer::AbstractWriter) = (writer.counter[] = 0)
reset_counter!(::Nothing) = nothing

##############################
# VTK Writer
##############################

"""
    VTKWriter(filename, mesh, fields)

"""
struct VTKWriter{N, T, C} <: AbstractWriter
    filename::String
    collectionfile::WriteVTK.CollectionFile

    coordinates::Vector{SVector{N, Float64}}
    shapes::Vector{T}
    dofprobes::Vector{PointProbe{N}}
    centerprobes::Vector{PointProbe{N}}

    vtkpoints::Vector{SVector{3, Float64}}
    vtkcells::Vector{MeshCell{VTKCellType, Vector{Int64}}}

    timestamps::Vector{Float64}
    counter::RefValue{Int}
    disabled::RefValue{Bool}

    cache::C
end

function VTKWriter(filename::String, pdes, fieldhandler;
        order=1, initcounter=0, disabled=false)

    collectionfile = paraview_collection(filename, append=false)

    fields = fieldhandler.fields
    field, = fields
    mesh = field.mesh
    domains = collect(eachdomaintag(mesh))
    coordinates = mesh.coordinates
    shapes = mesh.elements

    probesdofmap = DofMap(Lagrange(order), shapes; renumber=false)
    append!(probesdofmap.tags, domains)

    vtkpoints = pad3D.(get_dof_coordinates(mesh, probesdofmap))
    vtkcells = get_vtk_cells(mesh, probesdofmap)

    dofprobes = get_dof_pointprobes(mesh, probesdofmap)
    centerprobes = get_center_probes(mesh, domains)

    cache = OutputCache(pdes, fields)
    resize_cache1!(cache, length(dofprobes))
    resize_cache2!(cache, length(centerprobes))

    return VTKWriter{ndims(mesh), eltype(shapes), typeof(cache)}(
        filename,
        collectionfile,
        coordinates,
        shapes,
        dofprobes,
        centerprobes,
        vtkpoints,
        vtkcells,
        Float64[],
        Ref(initcounter),
        Ref(disabled),
        cache
    )
end

# TODO: Replace with Vector{DofProbe}?
# - No, because the discretization can be different between field and output? Would need new proxy and dofs.
function get_dof_pointprobes(mesh::AbstractMesh, dofmap::AbstractDofMap)
    probes = PointProbe{ndims(mesh)}[]
    resize!(probes, ndofs(dofmap))

    shapes = eachelement(mesh, dofmap.tags)

    return get_dof_pointprobes(dofmap.shapefunc, shapes, dofmap.dofs)
end

function get_dof_pointprobes(shapefunc::AbstractShapeFunctions, shapes, shapedofs)
    N = ndims(typeof(first(shapes)))
    n_dofs = maxdof(shapedofs)

    unassigned = trues(n_dofs)
    probes = Vector{PointProbe{N}}(undef, n_dofs)

    for (i, (shape, dofs)) in enumerate(zip(shapes, shapedofs))
        for (q, dof) in enumerate(dofs)
            if unassigned[dof]
                ξ = dof_refcoordinates(typeof(shape), shapefunc, q)
                probes[dof] = PointProbe(shape.id, ξ)
            end
        end
    end

    return probes
end

# TODO: Obsolete with ConstantShapeFunction?
function get_center_probes(mesh::AbstractMesh, domains=eachdomaintag(mesh))
    probes = PointProbe{ndims(mesh)}[]

    for idx in eachelementindex(mesh, domains)
        shape = mesh[idx]
        ξ = center_refcoordinates(typeof(shape))
        push!(probes, PointProbe(idx.shapeid, ξ))
    end

    return probes
end

function get_dof_coordinates(mesh::AbstractMesh, dofmap::AbstractDofMap)
    domains = dofmap.tags

    dofcoordinates = zeros(SVector{ndims(mesh), Float64}, ndofs(dofmap))
    proxy = Proxy{ndims(mesh)}(dofmap.shapefunc)

    for idx in eachelementindex(mesh, domains)
        shape = mesh[idx]
        dofs = dofmap[idx.shapeid]
        shapecoordinates = @view mesh.coordinates[shape.n]
        set_element_dof_evaluation!(proxy, typeof(shape), shapecoordinates)

        for (q, dof) in enumerate(dofs)
            dofcoordinates[dof] = proxy.x[q]
        end
    end

    return dofcoordinates
end

function get_vtk_cells(mesh::AbstractMesh, dofmap::AbstractDofMap)
    domains = dofmap.tags

    vtkcells = [MeshCell(vtkcelltype(typeof(mesh[idx]), dofmap.shapefunc), dofmap[idx.shapeid])
                for idx in eachelementindex(mesh, domains)]

    return vtkcells
end

pad3D(v::SVector{1}) = SVector{3, eltype(v)}(v[1], 0, 0)
pad3D(v::SVector{2}) = SVector{3, eltype(v)}(v[1], v[2], 0)
pad3D(v::SVector{3}) = SVector{3, eltype(v)}(v[1], v[2], v[3])

vtkcelltype(::Type{<:AbstractEdge})         = VTKCellTypes.VTK_LAGRANGE_CURVE
vtkcelltype(::Type{<:AbstractTriangle})     = VTKCellTypes.VTK_LAGRANGE_TRIANGLE
vtkcelltype(::Type{Quad4})                  = VTKCellTypes.VTK_QUAD
vtkcelltype(::Type{Quad8})                  = VTKCellTypes.VTK_QUADRATIC_QUAD
vtkcelltype(::Type{<:AbstractTetrahedron})  = VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON
vtkcelltype(::Type{Hex8})                   = VTKCellTypes.VTK_HEXAHEDRON
vtkcelltype(::Type{Hex20})                  = VTKCellTypes.VTK_QUADRATIC_HEXAHEDRON
vtkcelltype(::Type{Wedge6})                 = VTKCellTypes.VTK_WEDGE
vtkcelltype(::Type{Wedge15})                = VTKCellTypes.VTK_QUADRATIC_WEDGE

vtkcelltype(::Type{<:AbstractEdge},          ::Lagrange) = VTKCellTypes.VTK_LAGRANGE_CURVE
vtkcelltype(::Type{<:AbstractTriangle},      ::Lagrange) = VTKCellTypes.VTK_LAGRANGE_TRIANGLE
vtkcelltype(::Type{<:AbstractQuadrilateral}, ::Lagrange{1})     = VTKCellTypes.VTK_QUAD
vtkcelltype(::Type{<:AbstractQuadrilateral}, ::Lagrange{2})     = VTKCellTypes.VTK_QUADRATIC_QUAD
vtkcelltype(::Type{<:AbstractTetrahedron},   ::Lagrange) = VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON
vtkcelltype(::Type{<:AbstractHexahedron},    ::Lagrange{1})     = VTKCellTypes.VTK_HEXAHEDRON
vtkcelltype(::Type{<:AbstractHexahedron},    ::Lagrange{2})     = VTKCellTypes.VTK_QUADRATIC_HEXAHEDRON
vtkcelltype(::Type{<:AbstractWedge},         ::Lagrange{1})     = VTKCellTypes.VTK_WEDGE
vtkcelltype(::Type{<:AbstractWedge},         ::Lagrange{2})     = VTKCellTypes.VTK_QUADRATIC_WEDGE

# vtkcelltype_native(::Type{Edge2})   = VTKCellTypes.VTK_LINE
# vtkcelltype_native(::Type{Edge3})   = VTKCellTypes.VTK_QUADRATIC_EDGE
# vtkcelltype_native(::Type{Tri3})    = VTKCellTypes.VTK_TRIANGLE
# vtkcelltype_native(::Type{Tri6})    = VTKCellTypes.VTK_QUADRATIC_TRIANGLE
# vtkcelltype_native(::Type{Quad4})   = VTKCellTypes.VTK_QUAD
# vtkcelltype_native(::Type{Quad8})   = VTKCellTypes.VTK_QUADRATIC_QUAD
# vtkcelltype_native(::Type{Tet4})    = VTKCellTypes.VTK_TETRA
# vtkcelltype_native(::Type{Tet10})   = VTKCellTypes.VTK_QUADRATIC_TETRA
# vtkcelltype_native(::Type{Hex8})    = VTKCellTypes.VTK_HEXAHEDRON
# vtkcelltype_native(::Type{Hex20})   = VTKCellTypes.VTK_QUADRATIC_HEXAHEDRON
# vtkcelltype_native(::Type{Wedge6})  = VTKCellTypes.VTK_WEDGE
# vtkcelltype_native(::Type{Wedge15}) = VTKCellTypes.VTK_QUADRATIC_WEDGE

##############################
# Write
##############################

function write!(writer::VTKWriter, pdes, fieldhandler, t=0.0)
    writer.disabled[] && return

    pvd = writer.collectionfile

    timestamp = t == Inf ? 1 : t

    mkpath(dirname(writer.filename))

    pvd[timestamp] = writevtk(writer, pdes, fieldhandler)
    WriteVTK.save_file(pvd.xdoc, pvd.path) # TODO

    push!(writer.timestamps, t)
    writer.counter[] += 1

    return
end

# ------------------------------------------------------------

function writevtk(writer::VTKWriter, pdes, fieldhandler)
    filename = string(splitext(writer.filename)[1], "_", lpad(writer.counter[], 5, "0"))
    vtk = vtk_grid(filename, writer.vtkpoints, writer.vtkcells)

    appendfields!(vtk, writer, pdes, fieldhandler)

    return vtk
end

# ------------------------------------------------------------

function appendfields!(vtk, writer::VTKWriter, pdes, fieldhandler, mat=UndefinedMaterial())
    cache = writer.cache
    fields = fieldhandler.fields

    # VTK Point data
    @threads for i in eachindex(writer.dofprobes)
        probe = writer.dofprobes[i]
        setproxy!(fields, probe)
        append_nodal_fields!(cache.cache1, pdes, fields, mat, i)
    end

    # VTK Cell data
    @threads for i in eachindex(writer.centerprobes)
        probe = writer.centerprobes[i]
        setproxy!(fields, probe)
        append_elemental_fields!(cache.cache2, pdes, fields, mat, i)
    end

    # Looping is type instable
    for (name, field) in zip(cache.names1, cache.cache1)
        vtk[name] = field
    end

    for (name, field) in zip(cache.names2, cache.cache2)
        vtk[name] = field
    end
end

function _appendfields!(cache, pde::AbstractPDE, args...; kwargs...)
    _appendfields!(cache, (pde,), args...; kwargs...)
end

function append_nodal_fields!(cache::Tuple, pdes, fields, material, i)
    pde1, = pdes
    proxy = pde1.field.proxies[threadid()]
    x = proxy.x[1]

    vars = outputvars(pdes, fields, x, material)
    _set_each_field_index(cache, vars, i)
    return
end

function append_elemental_fields!(cache::Tuple, pdes, fields, material, i)
    pde1, = pdes
    proxy = pde1.field.proxies[threadid()]
    x = proxy.x[1]

    vars = outputvars_element(pdes, fields, x, material)
    _set_each_field_index(cache, vars, i)
    return
end

#############################
# Output Cache
#############################

struct OutputCache{T1, T2}
    cache1::T1
    names1::Vector{String}

    cache2::T2
    names2::Vector{String}
end

function OutputCache(pdes, fields)
    nodalvars = outputvars(_mock_output_args(pdes, fields)...)
    nodalvarnames = outputnames(pdes)

    cache1 = Tuple(Vector{typeof(field)}() for field in nodalvars)
    names1 = collect(getfield(nodalvarnames, key) for key in keys(nodalvars))

    elementvars = outputvars_element(_mock_output_args(pdes, fields)...)
    elementvarnames = outputnames_element(pdes)

    cache2 = Tuple(Vector{typeof(field)}() for field in elementvars)
    names2 = collect(getfield(elementvarnames, key) for key in keys(elementvars))

    return OutputCache(cache1, names1, cache2, names2)
end

# Looping over heterogeneoous Tuple is not type stable
function resize_cache1!(oc::OutputCache, n)
    for field in oc.cache1
        resize!(field, n)
    end
    return oc
end

function resize_cache2!(oc::OutputCache, n)
    for field in oc.cache2
        resize!(field, n)
    end
    return oc
end

function _mock_output_args(pdes, fields)
    testfield, = fields
    shape = first(eachshape(testfield))

    setproxy_point_evaluation!(fields, ElementIndex(shape.id), center_refcoordinates(typeof(shape)))

    x = first(testfield.mesh.coordinates)
    t = 0.0
    material = UndefinedMaterial()

    return pdes, fields, x, t, material
end

############################
# Default Output Variables
############################

outputvars(pde::AbstractPDE, args...) = outputvars((pde,), args...)

function outputvars(pdes::Tuple{Vararg{AbstractPDE, 1}}, args...)
    pde1, = pdes
    return (
        u1 = interpolate(pde1.field, 1),
    )
end

function outputvars(pdes::Tuple{Vararg{AbstractPDE, 2}}, args...)
    pde1, pde2 = pdes
    return (
        u1 = interpolate(pde1.field, 1),
        u2 = interpolate(pde2.field, 1),
    )
end

function outputvars(pdes::Tuple{Vararg{AbstractPDE, 3}}, args...)
    pde1, pde2, pde3 = pdes
    return (
        u1 = interpolate(pde1.field, 1),
        u2 = interpolate(pde2.field, 1),
        u3 = interpolate(pde3.field, 1),
    )
end

function outputvars(pdes::Tuple{Vararg{AbstractPDE, 4}}, args...)
    pde1, pde2, pde3, pde4 = pdes
    return (
        u1 = interpolate(pde1.field, 1),
        u2 = interpolate(pde2.field, 1),
        u3 = interpolate(pde3.field, 1),
        u4 = interpolate(pde4.field, 1),
    )
end

function outputvars(pdes::Tuple{Vararg{AbstractPDE, 5}}, args...)
    pde1, pde2, pde3, pde4, pde5, = pdes
    return (
        u1 = interpolate(pde1.field, 1),
        u2 = interpolate(pde2.field, 1),
        u3 = interpolate(pde3.field, 1),
        u4 = interpolate(pde4.field, 1),
        u5 = interpolate(pde5.field, 1),
    )
end

function outputvars(pdes::Tuple{Vararg{AbstractPDE, 6}}, args...)
    pde1, pde2, pde3, pde4, pde5, pde6 = pdes
    return (
        u1 = interpolate(pde1.field, 1),
        u2 = interpolate(pde2.field, 1),
        u3 = interpolate(pde3.field, 1),
        u4 = interpolate(pde4.field, 1),
        u5 = interpolate(pde5.field, 1),
        u6 = interpolate(pde6.field, 1),
    )
end

function outputvars(pdes::Tuple{Vararg{AbstractPDE, 7}}, args...)
    pde1, pde2, pde3, pde4, pde5, pde6, pde7 = pdes
    return (
        u1 = interpolate(pde1.field, 1),
        u2 = interpolate(pde2.field, 1),
        u3 = interpolate(pde3.field, 1),
        u4 = interpolate(pde4.field, 1),
        u5 = interpolate(pde5.field, 1),
        u6 = interpolate(pde6.field, 1),
        u7 = interpolate(pde7.field, 1),
    )
end

outputnames(pde::AbstractPDE, args...) = outputnames((pde,), args...)

function outputnames(pdes::Tuple{Vararg{AbstractPDE}})
    return NamedTuple(Symbol(:u, i) => pde.field.outputname for (i, pde) in enumerate(pdes))
end

# ------------------------------------------------------

outputvars_element(pde::AbstractPDE, args...) = outputvars_element((pde,), args...)

function outputvars_element(pdes::Tuple{Vararg{AbstractPDE, 1}}, args...)
    pde1, = pdes
    return (
        ∇u1 = interpolate_grad(pde1.field, 1),
    )
end

function outputvars_element(pdes::Tuple{Vararg{AbstractPDE, 2}}, args...)
    pde1, pde2 = pdes
    return (
        ∇u1 = interpolate_grad(pde1.field, 1),
        ∇u2 = interpolate_grad(pde2.field, 1),
    )
end

function outputvars_element(pdes::Tuple{Vararg{AbstractPDE, 3}}, args...)
    pde1, pde2, pde3 = pdes
    return (
        ∇u1 = interpolate_grad(pde1.field, 1),
        ∇u2 = interpolate_grad(pde2.field, 1),
        ∇u3 = interpolate_grad(pde3.field, 1),
    )
end

function outputvars_element(pdes::Tuple{Vararg{AbstractPDE, 4}}, args...)
    pde1, pde2, pde3, pde4 = pdes
    return (
        ∇u1 = interpolate_grad(pde1.field, 1),
        ∇u2 = interpolate_grad(pde2.field, 1),
        ∇u3 = interpolate_grad(pde3.field, 1),
        ∇u4 = interpolate_grad(pde4.field, 1),
    )
end

function outputvars_element(pdes::Tuple{Vararg{AbstractPDE, 5}}, args...)
    pde1, pde2, pde3, pde4, pde5, = pdes
    return (
        ∇u1 = interpolate_grad(pde1.field, 1),
        ∇u2 = interpolate_grad(pde2.field, 1),
        ∇u3 = interpolate_grad(pde3.field, 1),
        ∇u4 = interpolate_grad(pde4.field, 1),
        ∇u5 = interpolate_grad(pde5.field, 1),
    )
end

function outputvars_element(pdes::Tuple{Vararg{AbstractPDE, 6}}, args...)
    pde1, pde2, pde3, pde4, pde5, pde6 = pdes
    return (
        ∇u1 = interpolate_grad(pde1.field, 1),
        ∇u2 = interpolate_grad(pde2.field, 1),
        ∇u3 = interpolate_grad(pde3.field, 1),
        ∇u4 = interpolate_grad(pde4.field, 1),
        ∇u5 = interpolate_grad(pde5.field, 1),
        ∇u6 = interpolate_grad(pde6.field, 1),
    )
end

function outputvars_element(pdes::Tuple{Vararg{AbstractPDE, 7}}, args...)
    pde1, pde2, pde3, pde4, pde5, pde6, pde7 = pdes
    return (
        ∇u1 = interpolate_grad(pde1.field, 1),
        ∇u2 = interpolate_grad(pde2.field, 1),
        ∇u3 = interpolate_grad(pde3.field, 1),
        ∇u4 = interpolate_grad(pde4.field, 1),
        ∇u5 = interpolate_grad(pde5.field, 1),
        ∇u6 = interpolate_grad(pde6.field, 1),
        ∇u7 = interpolate_grad(pde7.field, 1),
    )
end

outputnames_element(pde::AbstractPDE, args...) = outputnames_element((pde,), args...)

function outputnames_element(pdes::Tuple{Vararg{AbstractPDE}})
    return NamedTuple(Symbol(:∇u, i) => ("∇" * pde.field.outputname * " [1/m]") for (i, pde) in enumerate(pdes))
end

#############################
# Misc
#############################

function _set_each_field_index(cache::NTuple{1, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 1
    c1, = cache
    v1, = vars
    c1[i] = v1
    return cache
end

function _set_each_field_index(cache::NTuple{2, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 2
    c1, c2 = cache
    v1, v2 = vars
    c1[i] = v1
    c2[i] = v2
    return cache
end

function _set_each_field_index(cache::NTuple{3, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 3
    c1, c2, c3 = cache
    v1, v2, v3 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    return cache
end

function _set_each_field_index(cache::NTuple{4, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 4
    c1, c2, c3, c4 = cache
    v1, v2, v3, v4 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    return cache
end

function _set_each_field_index(cache::NTuple{5, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 5
    c1, c2, c3, c4, c5 = cache
    v1, v2, v3, v4, v5 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    return cache
end

function _set_each_field_index(cache::NTuple{6, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 6
    c1, c2, c3, c4, c5, c6 = cache
    v1, v2, v3, v4, v5, v6 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    return cache
end

function _set_each_field_index(cache::NTuple{7, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 7
    c1, c2, c3, c4, c5, c6, c7 = cache
    v1, v2, v3, v4, v5, v6, v7 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    return cache
end

function _set_each_field_index(cache::NTuple{8, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 8
    c1, c2, c3, c4, c5, c6, c7, c8 = cache
    v1, v2, v3, v4, v5, v6, v7, v8 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    return cache
end

function _set_each_field_index(cache::NTuple{9, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 9
    c1, c2, c3, c4, c5, c6, c7, c8, c9 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    return cache
end

function _set_each_field_index(cache::NTuple{10, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 10
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9, v10 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    c10[i] = v10
    return cache
end

function _set_each_field_index(cache::NTuple{11, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 11
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    c10[i] = v10
    c11[i] = v11
    return cache
end

function _set_each_field_index(cache::NTuple{12, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 12
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    c10[i] = v10
    c11[i] = v11
    c12[i] = v12
    return cache
end

function _set_each_field_index(cache::NTuple{13, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 13
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    c10[i] = v10
    c11[i] = v11
    c12[i] = v12
    c13[i] = v13
    return cache
end

function _set_each_field_index(cache::NTuple{14, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 14
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    c10[i] = v10
    c11[i] = v11
    c12[i] = v12
    c13[i] = v13
    c14[i] = v14
    return cache
end

function _set_each_field_index(cache::NTuple{15, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 15
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    c10[i] = v10
    c11[i] = v11
    c12[i] = v12
    c13[i] = v13
    c14[i] = v14
    c15[i] = v15
    return cache
end

function _set_each_field_index(cache::NTuple{16, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 16
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    c10[i] = v10
    c11[i] = v11
    c12[i] = v12
    c13[i] = v13
    c14[i] = v14
    c15[i] = v15
    c16[i] = v16
    return cache
end

function _set_each_field_index(cache::NTuple{17, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 17
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    c10[i] = v10
    c11[i] = v11
    c12[i] = v12
    c13[i] = v13
    c14[i] = v14
    c15[i] = v15
    c16[i] = v16
    c17[i] = v17
    return cache
end

function _set_each_field_index(cache::NTuple{18, AbstractArray}, vars::Union{Tuple, NamedTuple}, i)
    @assert length(vars) == 18
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18 = cache
    v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18 = vars
    c1[i] = v1
    c2[i] = v2
    c3[i] = v3
    c4[i] = v4
    c5[i] = v5
    c6[i] = v6
    c7[i] = v7
    c8[i] = v8
    c9[i] = v9
    c10[i] = v10
    c11[i] = v11
    c12[i] = v12
    c13[i] = v13
    c14[i] = v14
    c15[i] = v15
    c16[i] = v16
    c17[i] = v17
    c18[i] = v18
    return cache
end
