export AbstractProbe
export PointProbe
export DofProbe
export getshape
export setproxy!
export locate
export to_refcoordinates
export PointProbes
export LineProbe
# export MeshProbe

abstract type AbstractProbe end

#############################
# Dof Probe
#############################

# Obsolete?
# struct DofProbe <: AbstractProbe
#     shapeid::Int
#     localdof::Int
# end

#############################
# Point Probe
#############################

"""
"""
struct PointProbe{N} <: AbstractProbe
    shapeid::Int
    refcoordinates::SVector{N, Float64}
end

Base.ndims(::Union{Type{<:PointProbe{N}}, PointProbe{N}}) where {N} = N

function PointProbe(shapeid::Integer, ξ::AbstractVector{<:Number})
    PointProbe(shapeid, SVector{length(ξ), Float64}(ξ))
end

"""
"""
function PointProbe{N}(x::AbstractVector{<:Number}, mesh::AbstractMesh,
        frame::AbstractFrame=SpatialFrame();
        inside_tol=1e-8, nonlinear_tol=1e-8) where {N}

    coordinates = getcoordinates(mesh, frame)
    shapes = getshapes(mesh, Val(N))
    e, ξ = locate(x, coordinates, shapes; error_on_fail=true, inside_tol, nonlinear_tol)
    return PointProbe(e, ξ)
end

function PointProbe(x::AbstractVector{<:Number}, mesh::AbstractMesh{N}, args...; kwargs...) where {N}
    @assert length(x) == N
    PointProbe{N}(x, mesh, args...; kwargs...)
end

"""
    PointProbes(points, mesh, frame=SpatialFrame) -> Vector{<:PointProbe}

TODO

If `error_on_fail=false` and `skip_failed=true`, a probe is skipped if the point could
not be located in the mesh withing the specified tolerances. Therefore, the returned vector
of probes may have less elements than the input vector of points.

If `error_on_fail=false` and `skip_failed=false`, a probe is created with `shapeid=0` and
`refcoordinates=NaN` if the point could not be located in the mesh withing the specified
tolerances. The returned vector of probes will have the same length as the input vector of
points.
"""
function PointProbes(points, mesh::AbstractMesh{N},
        frame::AbstractFrame=SpatialFrame();
        error_on_fail=true, skip_failed=true,
        inside_tol=1e-8, nonlinear_tol=1e-8) where {N}

    @assert length(first(points)) == N

    coordinates = getcoordinates(mesh, frame)
    shapes = getshapes(mesh, Val(N))
    es, ξs = locate(points, coordinates, shapes; error_on_fail, inside_tol, nonlinear_tol)

    if skip_failed
        probes = PointProbe{length(eltype(ξs))}[]
        for (e, ξ) in zip(es, ξs)
            e == 0 && continue
            push!(probes, PointProbe(e, ξ))
        end
    else
        probes = [PointProbe(e, ξ) for (e, ξ) in zip(es, ξs)]
    end

    return probes
end

@inline getshape(m::AbstractMesh, p::PointProbe{N}) where {N} = getshapes(m, Val(N))[p.shapeid]

"""
"""
function getdofpointprobes(field::AbstractField)
    shapes = eachshape(field)
    return _getdofpointprobes(field.shapefunc, shapes, field.dofs)
end

function getdofpointprobes(shapes, shapefunc::AbstractShapeFunctions; renumber=false)
    dofmap = DofMap(shapefunc, shapes; renumber)
    return getdofpointprobes(shapes, dofmap)
end

function getdofpointprobes(shapes, dofmap::AbstractDofMap)
    maxid = maximum(shape -> shape.id, shapes)
    dofs = Vector{Vector{Int}}(undef, maxid)
    for (shape, _dofs) in zip(shapes, dofmap.dofs)
        dofs[shape.id] = _dofs
    end
    return _getdofpointprobes(dofmap.shapefunc, shapes, dofs)
end

function _getdofpointprobes(shapefunc::AbstractShapeFunctions, shapes, dofs)
    N = ndims(typeof(first(shapes)))
    n_dofs = maxdof(dofs)

    unassigned = trues(n_dofs)
    probes = Vector{PointProbe{N}}(undef, n_dofs)

    for (i, shape) in enumerate(shapes)
        for (q, dof) in enumerate(dofs[shape.id])
            if unassigned[dof]
                ξ = dof_refcoordinates(typeof(shape), shapefunc, q)
                probes[dof] = PointProbe(shape.id, ξ)
                unassigned[dof] = false
            end
        end
    end
    return probes
end

function to_globalcoordinates(probe::PointProbe, mesh::AbstractMesh,
        frame::AbstractFrame=SpatialFrame())

    coordinates = getcoordinates(mesh, frame)
    shape = getshape(mesh, probe)
    shapecoordinates = @view coordinates[shape.n]
    return to_globalcoordinates(probe.refcoordinates, typeof(shape), shapecoordinates)
end

function to_globalcoordinates(probes::Vector{<:PointProbe}, mesh::AbstractMesh, args...)
    x = similar(mesh.coordinates, length(probes))
    to_globalcoordinates!(x, probes, mesh, args...)
end

function to_globalcoordinates!(target, probes::Vector{<:PointProbe}, mesh::AbstractMesh{N},
        frame::AbstractFrame=SpatialFrame()) where {N}

    coordinates = getcoordinates(mesh, frame)
    shapes = getshapes(mesh, Val(N))
    to_globalcoordinates!(target, probes, coordinates, shapes)
end

function to_globalcoordinates(probes::Vector{<:PointProbe}, coordinates, shapes::Vector{<:AbstractShape})
    x = similar(coordinates, length(probes))
    to_globalcoordinates!(x, probes, coordinates, shapes)
end

function to_globalcoordinates!(target, probes::Vector{<:PointProbe}, coordinates, shapes::Vector{<:AbstractShape})
    @assert length(target) == length(probes)
    for (i, probe) in enumerate(probes)
        shapetype = typeof(shapes[probe.shapeid])
        shapecoordinates = @view coordinates[shapes[probe.shapeid].n]
        target[i] = to_globalcoordinates(probe.refcoordinates, shapetype, shapecoordinates)
    end
    return target
end

function to_globalcoordinates_3D(probes::Vector{<:PointProbe}, mesh::AbstractMesh, args...)
    x = Vector{SVector{3, Float64}}(undef, length(probes))
    to_globalcoordinates!(x, probes, mesh, args...)
end

function to_globalcoordinates_3D!(target, probes::Vector{<:PointProbe}, mesh::AbstractMesh{N},
        frame::AbstractFrame=SpatialFrame()) where {N}

    coordinates = getcoordinates(mesh, frame)
    shapes = getshapes(mesh, Val(N))
    to_globalcoordinates_3D!(target, probes, coordinates, shapes)
end

function to_globalcoordinates_3D(probes::Vector{<:PointProbe}, coordinates, shapes::Vector{<:AbstractShape})
    x = Vector{SVector{3, Float64}}(undef, length(probes))
    to_globalcoordinates!(x, probes, coordinates, shapes)
end

function to_globalcoordinates_3D!(target, probes::Vector{<:PointProbe}, coordinates, shapes::Vector{<:AbstractShape})
    @assert length(target) == length(probes)
    for (i, probe) in enumerate(probes)
        shapetype = typeof(shapes[probe.shapeid])
        shapecoordinates = @view coordinates[shapes[probe.shapeid].n]
        target[i] = pad3D(to_globalcoordinates(probe.refcoordinates, shapetype, shapecoordinates))
    end
    return target
end

function setproxy!(field::AbstractDomainField, probe::PointProbe)
    setproxy_point_evaluation!(field, ElementIndex(probe.shapeid), probe.refcoordinates)
end

function setproxy!(field::AbstractBoundaryField, probe::PointProbe)
    setproxy_point_evaluation!(field, FacetIndex(probe.shapeid, 0, 0), probe.refcoordinates)
end

function setproxy!(notfield::Any, args...)
    return notfield
end

#############################
# Line Probe
#############################

struct LineProbe{N} <: AbstractProbe
    probes::Vector{PointProbe{N}}
end

"""
"""
function LineProbe(p1, p2, mesh::AbstractMesh{N},
        frame::AbstractFrame=SpatialFrame();
        n_points=500,
        error_on_fail=true, skip_failed=true,
        inside_tol=1e-8, nonlinear_tol=1e-8) where {N}

    @assert length(p1) == length(p2) == N

    points = [p1 + (p2 - p1) * i for i in range(0, 1, length=n_points)]
    points[1] = p1
    points[end] = p2

    probes = PointProbes(points, mesh, frame; error_on_fail, skip_failed, inside_tol, nonlinear_tol)
    return LineProbe{N}(probes)
end

function to_globalcoordinates(probe::LineProbe, mesh::AbstractMesh, args...)
    x = similar(mesh.coordinates, length(probe.probes))
    to_globalcoordinates!(x, probe, mesh, args...)
end

function to_globalcoordinates!(target, probe::LineProbe, mesh::AbstractMesh{N},
        frame::AbstractFrame=SpatialFrame()) where {N}

    coordinates = getcoordinates(mesh, frame)
    shapes = getshapes(mesh, Val(N))
    to_globalcoordinates!(target, probe.probes, coordinates, shapes)
end

function to_globalcoordinates_3D(probe::LineProbe, mesh::AbstractMesh, args...)
    x = Vector{SVector{3, Float64}}(undef, length(probe.probes))
    to_globalcoordinates!(x, probe, mesh, args...)
end

function to_globalcoordinates_3D!(target, probe::LineProbe, mesh::AbstractMesh{N},
        frame::AbstractFrame=SpatialFrame()) where {N}

    coordinates = getcoordinates(mesh, frame)
    shapes = getshapes(mesh, Val(N))
    to_globalcoordinates_3D!(target, probe.probes, coordinates, shapes)
end

#############################
# Locate Points
#############################

"""
"""
function locate(x::AbstractVector{<:Number}, coordinates, shapes;
        error_on_fail=true, inside_tol=1e-8, nonlinear_tol=1e-8, maxiter=100)

    @assert length(x) == length(first(coordinates))
    T = SVector{ndims(eltype(shapes)), Float64}

    vs = center_coordinates(coordinates, shapes) .- Ref(x)
    closest_index = findmin(v -> sum(v.^2), vs)[2]
    shapeneighbors = getneighbors(shapes, closest_index)

    shape_index = 0
    refcoordinates = NaN * zero(T)

    found = false
    for i in shapeneighbors
        shape = shapes[i]
        shapecoordinates = @view coordinates[shape.n]
        ξ = _to_refcoordinates(x, typeof(shape), shapecoordinates; tol=nonlinear_tol, maxiter)
        if isvalid_refcoordinates(ξ, typeof(shape), tol=inside_tol)
            shape_index = i
            refcoordinates = T(ξ)
            found = true
            break
        end
    end
    if !found && error_on_fail
        throw(ErrorException("could not locate point = $x in mesh"))
    end

    return shape_index, refcoordinates
end

"""

Returns index 0 and NaN-coordinates if `error_on_fail=false` and a point could not
be lcoated.
"""
function locate(points, coordinates, shapes;
        error_on_fail=true, inside_tol=1e-8, nonlinear_tol=1e-8, maxiter=100)

    @assert length(first(points)) == length(first(coordinates))
    T = SVector{ndims(eltype(shapes)), Float64}

    shapecenters = center_coordinates(coordinates, shapes)
    shapeneighbors = getneighbors(shapes)

    tree = NearestNeighbors.KDTree(shapecenters)
    nearest_shape, = NearestNeighbors.nn(tree, collect(points))
    @assert length(nearest_shape) == length(points)

    shape_indices = zeros(Int, length(points))
    refcoordinates = fill(NaN * zero(T), length(points))

    for (i, x) in enumerate(points)
        found = false
        for shape_index in shapeneighbors[nearest_shape[i]]
            shape = shapes[shape_index]
            shapecoordinates = @view coordinates[shape.n]
            ξ = _to_refcoordinates(x, typeof(shape), shapecoordinates; tol=nonlinear_tol, maxiter)
            if isvalid_refcoordinates(ξ, typeof(shape), tol=inside_tol)
                shape_indices[i] = shape_index
                refcoordinates[i] = ξ
                found = true
                break
            end
        end
        if !found && error_on_fail
            throw(ErrorException("could not locate points[$i] = $x in mesh"))
        end
    end
    return shape_indices, refcoordinates
end

"""
    to_refcoordinates(x, shapetype, shapecoordinates; tol=1e-8, maxiter=100)

Converts global coordinates `x` to reference coordinates `ξ` in `shapetype` with
`shapecoordinates`.

Uses Newton's method is case of nonlinear parametrization. In case of
`ndims(shapetype) < length(x)`, a Gauss-Newton method is used to solve
the least squares problem.

`maxiter`: Maximum number of Newton iterations.

`tol`: Convergece criteria (absolute tolerance ||x - x'(ξ)|| < tol)
"""
function to_refcoordinates(x, ::Type{S}, shapecoordinates; kwargs...) where {S <: AbstractShape}
    ξ = _to_refcoordinates(x, S, shapecoordinates; kwargs...)
    if !isvalid_refcoordinates(ξ, S)
        @warn "Could not locate $x inside the shape: Ref. coordinates ξ = $ξ"
    end
    return ξ
end

"""
    _to_refcoordinates(x, shapetype, shapecoordinates; tol=1e-8, maxiter=100)

Same as `to_refcoordinates`, but does not warn if the returned reference coordinate are outside the shape.
"""
function _to_refcoordinates(x, ::Type{S}, shapecoordinates; tol=1e-8, maxiter=100) where {S <: AbstractShape}
    @assert length(shapecoordinates) == ndofs(S)
    @assert length(x) == length(first(shapecoordinates))

    M = ndims(S)
    ξ = zero(SVector{M, Float64})
    tol2 = tol^2 * length(x)

    for k in 1:maxiter
        r = to_globalcoordinates(ξ, S, shapecoordinates) - x
        rnorm2 = sum(abs2, r)
        if rnorm2 < tol2
            break
        end
        drdξ = parameterization_jacobian(ξ, S, shapecoordinates)
        Δξ = _backslash_least_squares(drdξ, r)
        ξ -= Δξ
        if k == maxiter
            @warn "max. iterations to convert point $x to refcoordinates in $S reached"
        end
    end
    return ξ
end

function _to_refcoordinates(x, ::Type{Edge2}, shapecoordinates; kwargs...)
    @assert length(shapecoordinates) == 2
    @assert length(x) == length(first(shapecoordinates))

    p1, p2 = shapecoordinates
    J = SMatrix{length(x), 1, Float64}(p2 - p1)
    ξ = _backslash_least_squares(J, x .- p1)

    return ξ .- 0.5
end

function _to_refcoordinates(x, ::Type{Tri3}, shapecoordinates; kwargs...)
    @assert length(shapecoordinates) == 3
    @assert length(x) == length(first(shapecoordinates))

    p1, p2, p3 = shapecoordinates
    _p2 = p2 .- p1
    _p3 = p3 .- p1

    J = SMatrix{length(x), 2, Float64}(_p2..., _p3...)
    ξ = _backslash_least_squares(J, x .- p1)

    return ξ
end

function _to_refcoordinates(x, ::Type{Tet4}, shapecoordinates; kwargs...)
    @assert length(shapecoordinates) == 4
    @assert length(x) == length(first(shapecoordinates))

    p1, p2, p3, p4 = shapecoordinates
    _p2 = p2 .- p1
    _p3 = p3 .- p1
    _p4 = p4 .- p1

    J = SMatrix{length(x), 3, Float64}(_p2..., _p3..., _p4...)
    ξ = _backslash_least_squares(J, x .- p1)

    return ξ
end

# Dispatch to handwritten ordinary least squares solution for rectangular matrix A,
# because the backslash operator allocates when choosing the least squares solution
# even with SMatrix.
_backslash_least_squares(A::StaticMatrix, b) = (A' * A) \ (A' * b)
_backslash_least_squares(A::StaticMatrix{N, N}, b) where {N} = A \ b
_backslash_least_squares(A, b) = A \ b

"""
"""
function isvalid_refcoordinates end

function isvalid_refcoordinates(ξ, ::Type{S}; tol=1e-8) where {S <: Union{AbstractTriangle, AbstractTetrahedron}}
    @assert length(ξ) == ndims(S)
    return (sum(ξ) < 1.0 + tol * ndims(S)) && all(x -> -tol < x < 1.0 + tol, ξ)
end

function isvalid_refcoordinates(ξ, ::Type{S}; tol=1e-8) where {S <: Union{AbstractEdge, AbstractQuadrilateral, AbstractHexahedron}}
    @assert length(ξ) == ndims(S)
    return all(x -> -0.5 - tol < x < 0.5 + tol, ξ)
end

function isvalid_refcoordinates(ξ, ::Type{S}; tol=1e-8) where {S <: AbstractWedge}
    @assert length(ξ) == ndims(S)
    x, y, z = ξ
    return (-tol < x < 1.0 + tol) && (-tol < y < 1.0 + tol) && (x + y < 1.0 + tol * 2) && abs(z) < 0.5 + tol
end
