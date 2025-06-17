####################################
# Makie Convert Arguments Overload
####################################

function Makie.convert_arguments(::Type{<:Makie.Poly},
        coordinates::AbstractVector{<:AbstractVector},
        shapes::AbstractVector{<:AbstractShape2D})
    return (to_gb_mesh(coordinates, shapes), )
end

function Makie.convert_arguments(::Type{<:Makie.Mesh},
        coordinates::AbstractVector{<:AbstractVector},
        shapes::AbstractVector{<:AbstractShape2D})
    return (to_gb_mesh(coordinates, shapes), )
end

function Makie.convert_arguments(::Type{<:Makie.Wireframe},
        coordinates::AbstractVector{<:AbstractVector},
        shapes::AbstractVector{<:AbstractShape2D})
    return (to_gb_mesh(coordinates, shapes), )
end

function Makie.convert_arguments(::Type{<:Makie.Poly}, mesh::Union{AbstractMesh2D, AbstractMesh3D})
    return (to_gb_mesh(mesh), )
end

function Makie.convert_arguments(::Type{<:Makie.Mesh}, mesh::Union{AbstractMesh2D, AbstractMesh3D})
    return (to_gb_mesh(mesh), )
end

function Makie.convert_arguments(::Type{<:Makie.Wireframe}, mesh::Union{AbstractMesh2D, AbstractMesh3D})
    return (to_gb_mesh(mesh), )
end

function Makie.convert_arguments(T::Makie.PointBased, mesh::Union{AbstractMesh2D, AbstractMesh3D})
    return Makie.convert_arguments(T, mesh.coordinates)
end

function Makie.convert_arguments(T::Type{<:Makie.Spy}, soe::SystemOfEquations)
    return Makie.convert_arguments(T, soe.jacobian)
end

Makie.plottype(::AbstractVector{<:AbstractVector}, ::AbstractVector{<:AbstractShape2D}) = Makie.Poly

####################################
# Mesh Conversion to GeometryBasics
####################################

to_gbface(s::S) where {S <: AbstractShape2D} = gbfacetype(S)(firstorder(s).n...)
gbfacetype(::Type{<:AbstractTriangle}) = GeometryBasics.TriangleFace
gbfacetype(::Type{<:AbstractQuadrilateral}) = GeometryBasics.QuadFace

to_gbpoint(x::StaticVector{1}) = GeometryBasics.Point1f(x)
to_gbpoint(x::StaticVector{2}) = GeometryBasics.Point2f(x)
to_gbpoint(x::StaticVector{3}) = GeometryBasics.Point3f(x)

to_gbpoint_3f(x::StaticVector{1}) = GeometryBasics.Point3f(x[1], 0, 0)
to_gbpoint_3f(x::StaticVector{2}) = GeometryBasics.Point3f(x[1], x[2], 0)
to_gbpoint_3f(x::StaticVector{3}) = GeometryBasics.Point3f(x...)

"""

Converts Felder.jl `Mesh2D` or `Mesh3D` to a **first order**
GeometryBasics mesh. 3D meshes are returned as a boundary mesh.
"""
function to_gb_mesh(mesh::AbstractMesh2D, domaintags=eachdomaintag(mesh))
    elements = eachelement(mesh, domaintags)
    return to_gb_mesh(mesh.coordinates, elements)
end

function to_gb_mesh(mesh::AbstractMesh3D, domaintags=eachdomaintag(mesh))
    boundarytags = unique(boundarytag for domaintag in domaintags
        for boundarytag in mesh.domain2boundaries[domaintag])
    facets = eachfacet(mesh, boundarytags)
    return to_gb_mesh(mesh.coordinates, facets)
end

function to_gb_mesh(coordinates, shapes2D)
    points = to_gbpoint.(coordinates)
    faces = to_gbface.(shapes2D)
    if isconcretetype(eltype(faces))
        if eltype(points) <: GeometryBasics.Point{3}
            normals = GeometryBasics.normals(points, faces)
        else
            normals = fill(GeometryBasics.Vec3f(0, 0, 1), length(points))
        end
        return GeometryBasics.Mesh(GeometryBasics.meta(points; normals), faces)
    else
        return GeometryBasics.normal_mesh(points, faces)
    end
end

####################################
# Mesh Plot Recipe
####################################

"""
    meshlineplot(mesh1D)
    meshlineplot(coordinates, shapes1D)

## Attributes

### Specific to `MeshLinePlot`

- `xcomponent::Int = 1` sets the index of the spatial coordinate to be plotted
  on the x-axis (1, 2, or 3).

$(Base.Docs.doc(Makie.MakieCore.colormap_attributes!))

$(Base.Docs.doc(Makie.MakieCore.generic_plot_attributes!))
"""
Makie.@recipe(MeshLinePlot, coordinates, shapes1D) do scene
    attr = Makie.Attributes(;
        xcomponent=1,

        color=Makie.theme(scene, :linecolor),
        linewidth=Makie.theme(scene, :linewidth),

        marker=:circle,
        marker2=:circle,
        markersize=Makie.theme(scene, :markersize),
        markersize2=Makie.theme(scene, :markersize)[] * 0.75,

        linestyle=nothing,
        fxaa=false,
    )
    Makie.MakieCore.generic_plot_attributes!(attr)
    return Makie.MakieCore.colormap_attributes!(attr, Makie.theme(scene, :colormap))
end

function Makie.plot!(p::MeshLinePlot{<:Tuple{<:AbstractVector{<:AbstractVector}, <:AbstractVector{<:AbstractShape1D}}})
    coordinates = p.coordinates[]
    shapes = p.shapes1D[]
    xcomp = p.xcomponent[]

    vertices = SVector{2, Float64}[]
    edgenodes = SVector{2, Float64}[]
    for shape in shapes
        for i in 1:nvertices(typeof(shape))
            push!(vertices, SA[getindex(coordinates[shape.n[i]], xcomp), 0])
        end
        for i in nvertices(typeof(shape))+1:length(shape.n)
            push!(edgenodes, SA[getindex(coordinates[shape.n[i]], xcomp), 0])
        end
    end
    Makie.linesegments!(p, vertices;
        color=p.color,
        linestyle=p.linestyle,
        linewidth=p.linewidth,
        colormap=p.colormap,
        colorscale=p.colorscale,
        colorrange=p.colorrange,
        inspectable=p.inspectable
    )
    Makie.scatter!(p, vertices;
        color=p.color,
        marker=p.marker,
        markersize=p.markersize,
        inspectable=p.inspectable
    )
    Makie.scatter!(p, edgenodes;
        color=p.color,
        marker=p.marker2,
        markersize=p.markersize2,
        inspectable=p.inspectable
    )
    return p
end

function Makie.convert_arguments(::Type{<:MeshLinePlot}, mesh::AbstractMesh1D)
    coordinates = mesh.coordinates
    shapes = mesh.elements
    return (coordinates, shapes)
end

"""
    meshplot(mesh)
    meshplot(coordinates, shapes2D)

## Attributes

### Specific to `MeshPlot`

$(Base.Docs.doc(Makie.MakieCore.colormap_attributes!))

$(Base.Docs.doc(Makie.MakieCore.generic_plot_attributes!))
"""
Makie.@recipe(MeshPlot, mesh) do scene
    attr = Makie.Attributes(;
        color=:grey,
        alpha=1.0,
        strokecolor=:black,
        strokewidth=1,
        linestyle=nothing,
        shading=Makie.NoShading,
        fxaa=true,
    )
    return Makie.MakieCore.generic_plot_attributes!(attr)
end

function Makie.plot!(p::MeshPlot{<:Tuple{<:Union{AbstractMesh2D, AbstractMesh3D}}})
    Makie.poly!(p, p.mesh; p.attributes...)
    return p
end

function Makie.convert_arguments(::Type{<:MeshPlot},
        coordinates::AbstractVector{<:AbstractVector},
        shapes::AbstractVector{<:AbstractShape2D})
    return (to_gb_mesh(coordinates, shapes),)
end

Makie.plottype(::AbstractMesh1D) = MeshLinePlot
Makie.plottype(::AbstractMesh2D) = MeshPlot
Makie.plottype(::AbstractMesh3D) = MeshPlot

####################################
# Field Plot Recipe
####################################

"""
    fieldlineplot(field; kwargs...)

Line plot for one-dimensional fields.

## Attributes

### Specific to `FieldLinePlot`
- `component::Int = 0` sets the component of the field to be plotted if
  `field` is a `VectorField`. The default `component=0` plots the magnitude
  (2-norm) of the vector field. Ignored if `field` is a `ScalarField`.
- `xcomponent::Int = 1` sets the index of the spatial coordinate to be plotted
  on the x-axis (1, 2, or 3).
- `pointsperelement::Int = 10` sets the number of interpolation points per
  element (resolution).

$(Base.Docs.doc(Makie.MakieCore.generic_plot_attributes!))
"""
Makie.@recipe(FieldLinePlot, field) do scene
    attr = Makie.Attributes(;
        component=0,
        xcomponent=1,
        pointsperelement=10,
        markersize=0,
    )
    return attr
end

function Makie.plot!(p::FieldLinePlot{<:Tuple{<:AbstractField}})
    :label in keys(p.attributes) || (p.attributes[:label] = p.field[].outputname)
    xcomp = p.xcomponent[]
    comp = p.component[]

    x = Makie.Observable(Float64[])
    u = Makie.Observable(Float64[])

    ξs = collect(SA[ξ] for ξ in range(-0.5, 0.5, length=p.pointsperelement[]))
    S = eltype(collect(eachshape(p.field[])))

    function update_plot(field)
        empty!(x[])
        empty!(u[])
        initproxy_point_evaluation!(field, S, ξs)
        for index in eachshapeindex(field)
            setproxy_point_evaluation!(field, index)
            proxy = field.proxies[Threads.threadid()]
            for q in eachindex(proxy.x)
                push!(x[], getindex(proxy.x[q], xcomp))
                push!(u[], _compfunc(interpolate(field, q), comp))
            end
            push!(x[], NaN)
            push!(u[], NaN)
        end
        # Alternative implementation with PointProbe and single point evaluation
        # mesh = field.mesh
        # for shape in eachshape(field)
        #     for ξ in ξs
        #         probe = PointProbe(shape.id, ξ)
        #         setproxy!(field, probe)
        #         push!(x[], getindex(interpolate_coordinates(mesh, probe), xcomp))
        #         push!(u[], _compfunc(interpolate(field, 1), comp))
        #     end
        #     push!(x[], NaN)
        #     push!(u[], NaN)
        # end
        Makie.Observables.notify(x)
        Makie.Observables.notify(u)
    end

    Makie.Observables.onany(update_plot, p.field)
    update_plot(p.field[])

    Makie.scatterlines!(p, x, u; markersize=p.markersize)
    return p
end

_compfunc(u::Number, ::Integer) = u
_compfunc(u::AbstractVector, i::Integer) = (i == 0) ? norm(u) : getindex(u, i)

"""
    fieldplot(field; kwargs...)

For 2D and 3D fields. Fields are converted to first order. 3D fields are
only plotted on the boundaries.

Quads on Wedge boundaries are decomposed into triangles for compatibility
with GeometryBasics.Mesh.

## Attributes

### Specific to `FieldPlot`
- `component::Int = 0` sets the component of the field to be plotted if
  `field` is a `VectorField`. The default `component=0` plots the magnitude
  (2-norm) of the vector field. Ignored if `field` is a `ScalarField`.
- `frame::AbstractFrame = SpatialFrame()` sets the coordinate frame of the
   field. The default `SpatialFrame()` plots the field on the current mesh
   coordinates of the mesh , whereas `MaterialFrame()` plots the field on
   the initial mesh coordinates.
- `strokewidth::Real = 1` sets the width of the edges in the mesh.

$(Base.Docs.doc(Makie.MakieCore.colormap_attributes!))

$(Base.Docs.doc(Makie.MakieCore.generic_plot_attributes!))
"""
Makie.@recipe(FieldPlot, field) do scene
    attr = Makie.Attributes(;
        component=0,
        frame=SpatialFrame(),

        strokecolor=Makie.theme(scene, :patchstrokecolor),
        strokecolormap=Makie.theme(scene, :colormap),
        strokewidth=1,
        linestyle=nothing,

        shading=Makie.NoShading,
        fxaa=true,
    )
    Makie.MakieCore.generic_plot_attributes!(attr)
    return Makie.MakieCore.colormap_attributes!(attr, Makie.theme(scene, :colormap))
end

function Makie.plot!(p::FieldPlot{<:Tuple{<:AbstractField}})
    comp = p.component[]
    coordinates = getcoordinates(p.field[].mesh, p.frame[])
    shapes, x_indices, u_indices = Felder._convert2isoparametricfaces1(p.field[])
    @assert length(x_indices) == length(u_indices)

    u = Makie.Observable(zeros(length(u_indices)))
    gbmesh = Makie.Observable(to_gb_mesh(coordinates[x_indices], shapes))

    function update_plot(field)
        _u = u[]
        _gbmesh = gbmesh[]
        for i in eachindex(u_indices)
            _u[i] = _compfunc(field.u[u_indices[i]], comp)
            _gbmesh.position[i] = to_gbpoint(coordinates[x_indices[i]])
        end
        Makie.Observables.notify(u)
        Makie.Observables.notify(gbmesh)
    end

    Makie.Observables.onany(update_plot, p.field)
    update_plot(p.field[])

    delete!(p.attributes, :frame)
    delete!(p.attributes, :component)

    Makie.poly!(p, gbmesh; color=u, p.attributes...)
    return p
end

function Makie.plot!(p::FieldPlot{<:Tuple{<:AbstractField{<:Any, ConstantShapeFunction}}})
    comp = p.component[]
    coordinates = getcoordinates(p.field[].mesh, p.frame[])
    shapes, x_indices, u_indices = Felder._convert2faces1(p.field[])
    @assert length(shapes) == length(u_indices)

    coords = @view coordinates[x_indices]

    # Only u is updated as observable
    u = Makie.Observable(zeros(length(u_indices)))

    # This is a workaround because Makie.Poly does not plots polygons in 3D (z-coordinate
    # is dropped due to conversion to GeometryBasics.Point2f). Instead we plot the colored
    # polygons as a vector of meshes. However, Makie.wireframe still does not work
    # correctly (z-coordinate is dropped), but let's live with that for now (until
    # Makie.Poly and GeometryBasics support 3D polygons).
    # polygons = [GeometryBasics.Polygon(to_gbpoint.(coords[shape.n])) for shape in shapes]
    gbmeshes = Makie.Observable([GeometryBasics.normal_mesh(
        to_gbpoint.(coords[shape.n]),
        [to_gbface(typeof(shape)(1:length(shape.n)))]) for shape in shapes
    ])

    function update_plot(field)
        _u = u[]
        _gbmeshes = gbmeshes[]
        for (i, shape) in enumerate(shapes)
            _u[i] = _compfunc(field.u[u_indices[i]], comp)
            _gbmeshes[i].position .= to_gbpoint_3f.(coords[shape.n])
        end
        Makie.Observables.notify(u)
        Makie.Observables.notify(gbmeshes)
    end

    Makie.Observables.onany(update_plot, p.field)
    update_plot(p.field[])

    delete!(p.attributes, :component)
    delete!(p.attributes, :frame)

    delete!(p.attributes, :strokecolor)
    delete!(p.attributes, :strokecolormap)
    delete!(p.attributes, :linestyle)
    delete!(p.attributes, :strokewidth)

    # Makie.poly!(p, gbmeshes; color=u, p.attributes...)
    Makie.mesh!(p, gbmeshes; color=u, p.attributes...)

    return p
end

# function Makie.convert_arguments(::Type{<:FieldPlot},
#         coordinates::AbstractVector{<:AbstractVector},
#         shapes::AbstractVector{<:AbstractShape2D},
#         u::AbstractVector)
#     return (to_gb_mesh(coordinates, shapes), u)
# end

function Makie.plottype(field::AbstractField)
    if nrefdims(field) === 1
        return FieldLinePlot
    else
        return FieldPlot
    end
end

###########################################
# Field Surface Plot Recipe (for 2D Field)
###########################################

"""
    fieldsurface(field; kwargs...)
    fieldsurface(coordinates, shapes2D, u; kwargs...)

Surface plot 2D fields in 3D axes where the field variable `u` is
represented as the height of the surface.
Fields are converted to first order.

Not implemented yet for ConstantShapeFunction().

## Attributes

### Specific to `FieldPlot`
- `component::Int = 0` sets the component of the field to be plotted if
  `field` is a `VectorField`. The default `component=0` plots the magnitude
  (2-norm) of the vector field. Ignored if `field` is a `ScalarField`.
- `frame::AbstractFrame = SpatialFrame()` sets the coordinate frame of the
  field. The default `SpatialFrame()` plots the field on the current mesh
  coordinates of the mesh , whereas `MaterialFrame()` plots the field on
  the initial mesh coordinates.
- `strokewidth::Real = 1` sets the width of the edges in the mesh.

$(Base.Docs.doc(Makie.MakieCore.shading_attributes!))

$(Base.Docs.doc(Makie.MakieCore.colormap_attributes!))

$(Base.Docs.doc(Makie.MakieCore.generic_plot_attributes!))
"""
Makie.@recipe(FieldSurface, field) do scene
    attr = Makie.Attributes(;
        component=0,
        frame=SpatialFrame(),

        strokecolor=Makie.theme(scene, :patchstrokecolor),
        strokecolormap=Makie.theme(scene, :colormap),
        strokewidth=1,
        linestyle=nothing,

        shading=Makie.NoShading,
        fxaa=true,
    )
    Makie.MakieCore.generic_plot_attributes!(attr)
    return Makie.MakieCore.colormap_attributes!(attr, Makie.theme(scene, :colormap))
end

function Makie.plot!(p::FieldSurface{<:Tuple{<:AbstractField}})
    comp = p.component[]
    coordinates = getcoordinates(p.field[].mesh, p.frame[])
    shapes, x_indices, u_indices = Felder._convert2isoparametricfaces1(p.field[])
    @assert length(x_indices) == length(u_indices)

    u = zeros(length(u_indices))
    xyz = [SA[xy[1], xy[2], 0] for xy in @view coordinates[x_indices]]
    gbmesh = Makie.Observable(to_gb_mesh(xyz, shapes))

    function update_plot(field)
        for i in eachindex(u_indices)
            _u = _compfunc(field.u[u_indices[i]], comp)
            u[i] = _u

            x, y, = coordinates[x_indices[i]]
            xyz[i] = SA[x, y, _u]
        end
        gbmesh[] = to_gb_mesh(xyz, shapes)
    end

    Makie.Observables.onany(update_plot, p.field)
    update_plot(p.field[])

    delete!(p.attributes, :component)
    delete!(p.attributes, :frame)

    Makie.poly!(p, gbmesh; color=u, p.attributes...)
    return p
end

function Makie.plot!(::FieldSurface{<:Tuple{<:AbstractField{<:Any, ConstantShapeFunction}}})
    error("not implemented yet for ConstantShapeFunction()")
end
