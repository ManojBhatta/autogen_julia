export closest_point_brute
export closest_point
export signed_distance_brute
export signed_distance
export penetration_distance
export signed_projected_distance

##############################
# Closest Point
##############################

"""
    closest_point_brute(x, shapetype, shapecoordinates, n=20) -> x_min, ξ_min

Approximate the closest point on a shape to a given point `x` using a brute force method:
Generate a grid of points on/in the shape, check the distance to `x` for each point,
and return the closest point `x_min` and its reference coordinates `ξ_min`.

`n` defines how many grid points are generated in a single dimension of the shape.
The total number of evaluated grid points depends on the shape type:

| Shape         | Number of Evaluated Points | With n = 20   |
|---------------|----------------------------|---------------|
| Edge          | n                          | 20            |
| Triangle      | sum(1:n)                   | 210           |
| Quadrilateral | n * n                      | 400           |
| Tetrahedron   | sum(i -> sum(1:i), 1:n)    | 1540          |
| Hexahedron    | n * n * n                  | 8000          |
| Wedge         | sum(1:n) * n               | 2100          |

Works for all shapes in all dimensions.

See also `closest_point` for a higher precision method using optimization.

# Example:
```julia-repl
julia> S = Tri6
Tri6

julia> shapecoordinates = [
 SA[0.0, 0.0, 0.2],
 SA[1.0, 0.0, 0.1],
 SA[0.0, 1.0, 0.2],
 SA[0.5, 0.0, 0.2],
 SA[0.5, 0.5, 0.0],
 SA[0.0, 0.5, 0.1]];

julia> x = SA[0.1, 0.2, 0.5];

julia> x_min, ξ_min = closest_point_brute(x, S, shapecoordinates, 20);

julia> x_min
3-element SVector{3, Float64} with indices SOneTo(3):
 0.10526315789473686
 0.15789473684210525
 0.14847645429362882

julia> ξ_min
2-element SVector{2, Float64} with indices SOneTo(2):
 0.10526315789473684
 0.15789473684210525

julia> @btime Felder.closest_point_brute(\$x, Tri6, \$shapecoordinates);
 3.438 μs (0 allocations: 0 bytes)
```
"""
function closest_point_brute(x, ::Type{S}, shapecoordinates, n=20) where {S <: AbstractShape}
    dist2_min = Inf
    x_min = zero(similar_type(x))
    ξ_min = zero(similar_type(x, Size(ndims(S))))
    for ξ in refcoordinates_grid(S, n)
        x_ξ = to_globalcoordinates(ξ, S, shapecoordinates)
        dist2 = sum(abs2, x - x_ξ)
        if dist2 < dist2_min
            dist2_min = dist2
            x_min = x_ξ
            ξ_min = ξ
        end
    end
    return x_min, ξ_min
end

"""
    closest_point_brute(x, coordinates, shapes) -> x_min, ξ_min, shapeindex

Generalization of `closest_point_brute` for multiple `shapes`. Return the closest point `x_min`,
its reference coordinates `ξ_min`, and the index of the closest shape in `shapes`.
"""
function closest_point_brute(x, coordinates, shapes::AbstractArray{S}, n=20) where {S <: AbstractShape}
    mindist2 = Inf
    x_min = zero(x)
    ξ_min = zero(SVector{ndims(S)})
    shapeindex = 0
    for (i, shape) in enumerate(shapes)
        shapecoordinates = @view coordinates[shape.n]
        x_ξ, ξ = closest_point_brute(x, typeof(shape), shapecoordinates,n)
        dist2 = sum(abs2, x_ξ - x)
        if dist2 < mindist2
            mindist2 = dist2
            x_min = x_ξ
            ξ_min = ξ
            shapeindex = i
        end
    end
    shapeindex != 0 || error("No closest point found")
    return x_min, ξ_min, shapeindex
end

"""
    closest_point(x, shapetype, shapecoordinates; kwargs...) -> x_min, ξ_min

Find the closest point `x_min` on a shape to a given point `x` by optimization
(nonlinear-least squares: minimzation of distance with Gauss-Newton method).
Works for 1D shapes in 2D and 3D space, and for 2D shapes in 3D space.

Note: Linear shapes `Edge2` and `Tri3` are using analytical solutions.

`opt_tol`: Convergence criterion for the Gauss-Newton method (gradient norm)

`maxiter`: The maximum number of iterations for the Gauss-Newton method

Checks all shape vertices as initial guesses.

See also `closest_point_brute` for a brute force method without optimization.

# Example:
```julia-repl
julia> S = Tri6
Tri6

julia> shapecoordinates = [
 SA[0.0, 0.0, 0.2],
 SA[1.0, 0.0, 0.1],
 SA[0.0, 1.0, 0.2],
 SA[0.5, 0.0, 0.2],
 SA[0.5, 0.5, 0.0],
 SA[0.0, 0.5, 0.1]];

julia> x = SA[0.1, 0.2, 0.5];

julia> x_min, ξ_min = closest_point(x, S, shapecoordinates);

julia> x_min
3-element SVector{3, Float64} with indices SOneTo(3):
 0.09986598736261608
 0.1523703537805299
 0.15024385769211568

julia> ξ_min
2-element SVector{2, Float64} with indices SOneTo(2):
 0.09986598736261608
 0.1523703537805299

julia> @btime Felder.closest_point(\$x, Tri6, \$shapecoordinates);
  1.910 μs (0 allocations: 0 bytes)

julia> x = SA[10.0, 20.0, 0.0];

julia> @btime Felder.closest_point(\$x, Tri6, \$shapecoordinates);
  3.375 μs (0 allocations: 0 bytes)
```
"""
function closest_point(x, ::Type{S}, shapecoordinates; kwargs...) where {S <: AbstractShape1D}
    @assert length(x) > 1

    mindist2 = Inf
    x_min = zero(x)
    ξ_min = zero(SVector{ndims(S)})
    for i in 1:nvertices(S)
        # Check distance with vertex
        ξ_vertex = dof_refcoordinates(S, sftype(S)(), i)
        x_vertex = to_globalcoordinates(ξ_vertex, S, shapecoordinates)
        dist2 = sum(abs2, x_vertex - x)
        if dist2 < mindist2
            mindist2 = dist2
            x_min = x_vertex
            ξ_min = ξ_vertex
        end

        # Minimize distance along edge with vertex as initial guess
        ξ = _closest_refcoordinates(x, S, shapecoordinates, ξ_vertex; kwargs...)
        if isvalid_refcoordinates(ξ, S; tol=0.0)
            x_ξ = to_globalcoordinates(ξ, S, shapecoordinates)
            dist2 = sum(abs2, x_ξ - x)
            if dist2 < mindist2
                mindist2 = dist2
                x_min = x_ξ
                ξ_min = ξ
            end
        end
    end

    return x_min, ξ_min
end

function closest_point(x, ::Type{S}, shapecoordinates; kwargs...) where {S <: AbstractShape2D}
    # TODO: Refactor
    @assert length(x) == 3

    mindist2 = Inf
    x_min = zero(x)
    ξ_min = zero(SVector{ndims(S)})
    for i in 1:nvertices(S)
        # Check distance with vertex
        ξ_vertex = dof_refcoordinates(S, sftype(S)(), i)
        x_vertex = to_globalcoordinates(ξ_vertex, S, shapecoordinates)
        dist2 = sum(abs2, x_vertex - x)
        if dist2 < mindist2
            mindist2 = dist2
            x_min = x_vertex
            ξ_min = ξ_vertex
        end

        # Minimize distance on face with vertex as initial guess
        ξ = _closest_refcoordinates(x, S, shapecoordinates, ξ_vertex; kwargs...)
        if isvalid_refcoordinates(ξ, S; tol=0.0)
            x_ξ = to_globalcoordinates(ξ, S, shapecoordinates)
            dist2 = sum(abs2, x_ξ - x)
            if dist2 < mindist2
                mindist2 = dist2
                x_min = x_ξ
                ξ_min = ξ
            end
        end
    end

    # Minimize distance on face with center as initial guess
    ξ = _closest_refcoordinates(x, S, shapecoordinates, center_refcoordinates(S); kwargs...)
    if isvalid_refcoordinates(ξ, S; tol=0.0)
        x_ξ = to_globalcoordinates(ξ, S, shapecoordinates)
        dist2 = sum(abs2, x_ξ - x)
        if dist2 < mindist2
            mindist2 = dist2
            x_min = x_ξ
            ξ_min = ξ
        end
    end

    for i in 1:nedges(S)
        S_edge = localedgetype(S, i)
        edgecoordinates = @view shapecoordinates[localedgeindices(S, i)]

        for j in 1:nvertices(S_edge)
            ξ_vertex = dof_refcoordinates(S_edge, sftype(S_edge)(), j)
            # Minimize distance along edge with vertex as initial guess
            ξ = _closest_refcoordinates(x, S_edge, edgecoordinates, ξ_vertex; kwargs...)
            if isvalid_refcoordinates(ξ, S_edge; tol=0.0)
                x_ξ = to_globalcoordinates(ξ, S_edge, edgecoordinates)
                dist2 = sum(abs2, x_ξ - x)
                if dist2 < mindist2
                    mindist2 = dist2
                    x_min = x_ξ
                    ξ_min = edge2facecoordinates(S, i, ξ)
                end
            end
        end
    end

    return x_min, ξ_min
end

function closest_point(x, ::Type{Edge2}, shapecoordinates; kwargs...)
    @assert length(shapecoordinates) == 2
    p1, p2 = shapecoordinates
    x_min, λ = closest_point_edge_2D_3D(x, p1, p2)
    ξ_min = SA[λ - 0.5]
    return x_min, ξ_min
end

function closest_point(x, ::Type{Tri3}, shapecoordinates; kwargs...)
    @assert length(shapecoordinates) == 3
    p1, p2, p3 = shapecoordinates
    x_min, ξ1, ξ2 = closest_point_triangle_3D(x, p1, p2, p3)
    ξ_min = SA[ξ1, ξ2]
    return x_min, ξ_min
end

"""
    closest_point(x, coordinates, shapes) -> x_min, ξ_min, shapeindex

Generalization of `closest_point` for multiple `shapes`. Return the closest point `x_min`,
its reference coordinates `ξ_min`, and the index of the closest shape in `shapes`.

Returns the first closest shape if multiple shapes have the same distance.
"""
function closest_point(x, coordinates, shapes::AbstractArray{S}; kwargs...) where {S <: AbstractShape}
    mindist2 = Inf
    x_min = zero(x)
    ξ_min = zero(SVector{ndims(S)})
    shapeindex = 0
    for (i, shape) in enumerate(shapes)
        shapecoordinates = @view coordinates[shape.n]
        x_ξ, ξ = closest_point(x, typeof(shape), shapecoordinates; kwargs...)
        dist2 = sum(abs2, x_ξ - x)
        if dist2 < mindist2
            mindist2 = dist2
            x_min = x_ξ
            ξ_min = ξ
            shapeindex = i
        end
    end
    shapeindex != 0 || error("No closest point found")
    return x_min, ξ_min, shapeindex
end

##############################
# Signed Distance
##############################

"""
    signed_distance_brute(x, shapetype, shapecoordinates, n=6) -> signed_dist, normal

Approximate the signed distance from a point `x` to a shape using a brute force method
and return it with the `normal` vector at the closest point on the shape.

The signed distance is negative if the dot product `dot(x - x_min, normal)` is negative,
where `x_min` is the closest point on the shape to `x` and `normal` is the normal vector
at `x_min`. Otherwise, the signed distance positive.
The normal vector on facets is defined as pointing outwards the element or domain.

Works for 1D shapes in 2D and 3D space, and for 2D shapes in 3D space.

See `closest_point_brute` for more information on the brute force method.

Note: Pay attention to the discontinuity when the sign switches with nonzero distance.
This can even happen **on** the shape because of the discretization error of the
brute force evaluation points.

`n` defines how many grid points are generated in a single dimension of the shape.
The total number of evaluated grid points depends on the shape type:

| Shape         | Number of Evaluated Points | With n = 20   |
|---------------|----------------------------|---------------|
| Edge          | n                          | 20            |
| Triangle      | sum(1:n)                   | 210           |
| Quadrilateral | n * n                      | 400           |
| Tetrahedron   | sum(i -> sum(1:i), 1:n)    | 1540          |
| Hexahedron    | n * n * n                  | 8000          |
| Wedge         | sum(1:n) * n               | 2100          |
"""
function signed_distance_brute(x, ::Type{S}, shapecoordinates, n=20) where {S <: AbstractShape}
    x_min, ξ_min = closest_point_brute(x, S, shapecoordinates, n)
    d = x - x_min
    normal = parameterization_normal(ξ_min, S, shapecoordinates)
    normal /= norm(normal)
    signed_dist = dot(d, normal) >= 0 ? norm(d) : -norm(d)
    return signed_dist, normal
end

"""
    signed_distance_brute(x, coordinates, shapes) -> signed_dist, normal

Generalization of `signed_distance_brute` for multiple `shapes`.
"""
function signed_distance_brute(x, coordinates, shapes::AbstractArray{S}, n=20) where {S <: AbstractShape}
    mindist2 = Inf
    x_mean = zero(x)
    normal_mean = zero(x)
    counter = 0
    for shape in shapes
        shapecoordinates = @view coordinates[shape.n]
        x_ξ, ξ = closest_point_brute(x, typeof(shape), shapecoordinates, n)
        d = x - x_ξ
        dist2 = sum(abs2, d)
        if dist2 > mindist2
            continue
        end
        normal = parameterization_normal(ξ, typeof(shape), shapecoordinates)
        normal = normal / norm(normal)
        if dist2 == mindist2 # Same minimum -> Average normal
            normal_mean += normal
            x_mean += x_ξ
            counter += 1
        else # New minimum
            mindist2 = dist2
            normal_mean = normal
            x_mean = x_ξ
            counter = 1
        end
    end
    x_mean /= counter
    minnormal = normal_mean / norm(normal_mean)
    mindist = dot(x - x_mean, minnormal) < 0 ? -sqrt(mindist2) : sqrt(mindist2)
    return mindist, minnormal
end

"""
    signed_distance(x, shapetype, shapecoordinates; kwargs...) -> signed_dist, normal

Compute the signed distance from a point `x` to a shape using optimization
(nonlinear-least squares: minimzation of distance with Gauss-Newton method) and
return it with the `normal` vector at the closest point on the shape.

The signed distance is negative if the dot product `dot(x - x_min, normal)` is negative,
where `x_min` is the closest point on the shape to `x` and `normal` is the normal vector
at `x_min`. Otherwise, the signed distance positive.
The normal vector on facets is defined as pointing outwards the element or domain.

Works for 1D shapes in 2D and 3D space, and for 2D shapes in 3D space.

See `closest_point` for more information on the optimization method.

Note: Pay attention to the discontinuity when the sign switches with nonzero distance
(when moving **around** the shape).
"""
function signed_distance(x, ::Type{S}, shapecoordinates; kwargs...) where {S <: AbstractShape}
    x_min, ξ_min = closest_point(x, S, shapecoordinates; kwargs...)
    d = x - x_min
    normal = parameterization_normal(ξ_min, S, shapecoordinates)
    normal /= norm(normal)
    signed_dist = dot(d, normal) >= 0 ? norm(d) : -norm(d)
    return signed_dist, normal
end

"""
    signed_distance(x, coordinates, shapes) -> signed_dist, normal

Generalization of `signed_distance` for multiple `shapes`.
"""
function signed_distance(x, coordinates, shapes::AbstractArray{S}; kwargs...) where {S <: AbstractShape}
    mindist2 = Inf
    x_mean = zero(x)
    normal_mean = zero(x)
    counter = 0
    for shape in shapes
        shapecoordinates = @view coordinates[shape.n]
        x_ξ, ξ = closest_point(x, typeof(shape), shapecoordinates; kwargs...)
        d = x - x_ξ
        dist2 = sum(abs2, d)
        if dist2 > mindist2
            continue
        end
        normal = parameterization_normal(ξ, typeof(shape), shapecoordinates)
        normal = normal / norm(normal)
        if dist2 == mindist2 # Same minimum -> Average normal
            normal_mean += normal
            x_mean += x_ξ
            counter += 1
        else # New minimum
            mindist2 = dist2
            normal_mean = normal
            x_mean = x_ξ
            counter = 1
        end
    end
    x_mean /= counter
    minnormal = normal_mean / norm(normal_mean)
    mindist = dot(x - x_mean, minnormal) < 0 ? -sqrt(mindist2) : sqrt(mindist2)
    return mindist, minnormal
end

##############################
# Projected Distance
##############################

"""
    signed_projected_distance(x, dir, shapetype, shapecoordinates; kwargs...) -> dist

Compute the signed projected distance `dist` from a point `x` to the intersection
point on the shape along direction `dir`. `dist` is positive if `x`
is on the same side in which the shape normal at the intersection point points,
otherwise negative.

Return `NaN` if the ray `x + t * dir/norm(dir)` with parameter `T` does
not intersect the shape.

Works for 1D shapes in 2D, and for 2D shapes in 3D space.

Note: The projected distance is not necessarily the minimum distance to the shape,
and `x_proj` is not necessarily the closest point on the shape to `x`.
"""
function signed_projected_distance(x, dir, ::Type{S}, shapecoordinates; kwargs...) where {S <: AbstractShape}
    normdir = dir / norm(dir)
    _signed_projected_distance(x, normdir, S, shapecoordinates; kwargs...)
end

function _signed_projected_distance(x, normdir, ::Type{S}, shapecoordinates; kwargs...) where {S <: AbstractShape}
    ξ, t = projected_point(x, normdir, S, shapecoordinates; kwargs...)
    if isnan(t) || iszero(t)
        dist = t
    else
        n = parameterization_normal(ξ, S, shapecoordinates)
        normal = n / norm(n)
        a = dot(normal, sign(t) * normdir)
        if abs(a) < 0.1 # Reject tangent rays because sign is ambiguous
            dist = NaN
        else
            dist = a < 0 ? abs(t) : -abs(t)
        end
    end
    return dist
end

"""
    signed_projected_distance(x, dir, coordinates, shapes; kwargs...) -> dist, shapeindex

Generalization of `signed_projected_distance` for multiple `shapes`.

Return `dist=NaN` and `shapeindex=0` if `x + t * dir/norm(dir)` with parameter `t`
does not intersect with any shape.
"""
function signed_projected_distance(x, dir, coordinates, shapes::AbstractArray{S}; kwargs...) where {S <: AbstractShape}
    normdir = dir / norm(dir)
    _signed_projected_distance(x, normdir, coordinates, shapes; kwargs...)
end

function _signed_projected_distance(x, dir, coordinates, shapes::AbstractArray{S}; kwargs...) where {S <: AbstractShape}
    minabsdist = Inf
    signed_dist = NaN
    shapeindex = 0
    for (i, shape) in enumerate(shapes)
        shapecoordinates = @view coordinates[shape.n]
        dist = _signed_projected_distance(x, dir, typeof(shape), shapecoordinates; kwargs...)
        if abs(dist) < minabsdist
            minabsdist = abs(dist)
            signed_dist = dist
            shapeindex = i
        end
    end
    return signed_dist, shapeindex
end

##############################
# Projected Point
##############################

"""
    projected_point(x, dir, shapetype, shapecoordinates, ξ_init=shapecenter; tol=1e-8) -> ξ, t

Find the reference coordinates `ξ` of the point on a shape intersecting the
ray `x + t * dir` with scalar parameter `t` using Newton's method.

In other words, return `ξ` and `t` such that

    to_globalcoordinates(ξ, shapetype, shapecoordinates) ≈ x + t * dir

within a tolerance `tol` (L2-norm).

The initial guess for the reference coordinates `ξ` can be passed with `ξ_init`
(default: shape center). The initial guess for `t` is zero.

Return `ξ=NaN` and `t=NaN` if the Newton method does not converge within `maxiter` or
if the determined reference coordinates `ξ` are outside the shape (using
`isvalid_refcoordinates` with tolerance `inside_tol`).

Note: Tangent rays lead to bad convergence of the Newton solver. The default
`maxiter=10` is chosen to avoid long computation times and to "detect" such rays
simply by returning `NaN` early.

Works for 1D shapes in 2D, and for 2D shapes in 3D space.

Note: In theory there are multiple solutions possible for higher order geometry,
but only the one found by Newton's method with the above initial guesses is returned.
"""
function projected_point(x, dir, ::Type{S}, shapecoordinates, ξ_init=center_refcoordinates(S);
        tol=1e-8, inside_tol=1e-8, maxiter=10) where {S <: AbstractShape}

    @assert length(x) == length(dir)
    @assert length(x) == ndims(S) + 1
    @assert length(shapecoordinates) == ndofs(S)
    @assert length(x) == length(first(shapecoordinates))

    tol2 = tol^2 * ndims(S)
    t = 0.0
    ξ = deepcopy(ξ_init)
    for k in 1:maxiter
        r = to_globalcoordinates(ξ, S, shapecoordinates) - x - t * dir
        if sum(abs2, r) < tol2
            if !isvalid_refcoordinates(ξ, S, tol=inside_tol)
                t *= NaN
                ξ *= NaN
            end
            break
        end
        drdt = -dir
        drdξ = parameterization_jacobian(ξ, S, shapecoordinates)
        J = hcat(drdt, drdξ)
        Δt, Δξ... = J \ -r
        t += Δt
        ξ += Δξ
        if k == maxiter
            t *= NaN
            ξ *= NaN
        end
    end
    return ξ, t
end

##############################
# Minimum Distance
##############################

"""
    _closest_refcoordinates(x, shapetype, shapecoordinates, ξ_init=shapecenter; kwargs...) -> ξ

Find the reference coordinates `ξ` on a shape that minimize the distance to a given point `x`
using Newton's method with initial guess `ξ_init`. (unconstrained optimization)

Returns `ξ=NaN` if the Newton method does not converge within `maxiter` iterations or
if the Hessian matrix at `ξ` is not positive definite (using Sylvester's criterion).

Reference coordinates `ξ` are not constrained to the shape. Use `isvalid_refcoordinates`
to check the validity of the reference coordinates.

Convergence is reached if the L2-norm of the gradient of the squared distance function
is less than `tol`.

Works for 1D shapes in 2D and 3D, and for 2D shapes in 3D space.

Note: Saddle and inflection points in the distance function lead to bad convergence
of the Newton solver. The default `maxiter=10` is chosen to avoid long computation times
and to "detect" such points simply by returning `ξ=NaN` early.

Note: This algorithm convergeces to a local minimum, which is not necessarily the global minimum.
Use with different initial guesses `ξ_init` to find the global minimum on the shape for higher
order geometry. The minimum distance can also be to a vertex or to an edge of the shape.
"""
function _closest_refcoordinates(x, ::Type{S}, shapecoordinates, ξ_init=center_refcoordinates(S);
        tol=1e-8, maxiter=10) where {S <: AbstractShape}

    @assert length(shapecoordinates) == ndofs(S)
    @assert length(x) == length(first(shapecoordinates))

    tol2 = tol^2 * length(x)
    ξ = deepcopy(ξ_init)
    for k in 1:maxiter
        r = to_globalcoordinates(ξ, S, shapecoordinates) - x
        drdξ = parameterization_jacobian(ξ, S, shapecoordinates)
        # f = 0.5 * sum(abs2, r)
        ∇f = drdξ' * r
        ∇f_norm2 = sum(abs2, ∇f)
        ∇∇f = drdξ' * drdξ # Gauss-Newton
        ∇∇f += hesse_higher_order_term(r, ξ, S, shapecoordinates)
        if ∇f_norm2 < tol2
            if ∇∇f[1, 1] <= 0 || det(∇∇f) <= 0 # Sylvester's criterion
                ξ *= NaN
            end
            break
        end
        Δξ = ∇∇f \ -∇f # Newton direction
        ξ = ξ + Δξ
        if k == maxiter
            ξ *= NaN
        end
    end
    return ξ
end

function hesse_higher_order_term(r, ξ, ::Type{S}, shapecoordinates) where {S <: AbstractShape}
    @assert ndofs(S) == length(shapecoordinates)
    @assert ndims(S) == length(ξ)
    @assert length(first(shapecoordinates)) == length(r)
    M = ndims(S)
    H_term = zero(SMatrix{M, M})
    for (xi, ∇∇Ni) in zip(shapecoordinates, eachshapefunc_∇∇N(S, sftype(S)(), ξ))
        for (xij, rj) in zip(xi, r)
            H_term += xij * ∇∇Ni * rj
        end
    end
    return H_term
end
