export interpolate
export interpolate!
export interpolate_grad!
export interpolate_nodalmean!
export interpolate_boundary!
export interpolate_edge!  #TODO
export interpolate_vertex! # TODO
export eval_interpolant
export SpaceFunction
export TimeFunction
export SpaceTimeFunction

# TODO: Vector components
# TODO: Checks that the tags actually exists in the field

#################################
# Dispatch on function arguments
#################################

"""
    f = SpaceFunction(func::Function)

1-argument function of space coordinate `x` that can be evaluated as `f(x)`.

Can be used in functions like `interpolate!` and `integrate` to dispatch
on the intended input arguments.

# Example:
```julia-repl
julia> f = SpaceFunction(x -> sin(x[1]));

julia> f(π/2)
1.0

julia> f([1.5π, 0, 0])
-1.0

julia> g = SpaceFunction(x -> 2 * x);

julia> g([1, 2, 3])
3-element Vector{Int64}:
 2
 4
 6
```
"""
struct SpaceFunction{F}
    f::F
end

(f::SpaceFunction{F})(x) where {F} = f.f(x)
(f::SpaceFunction{F})(x) where {F <: Number} = f.f
(f::SpaceFunction{F})(x) where {F <: AbstractArray} = f.f
(f::SpaceFunction{F})(x) where {F <: RefValue} = f.f[]

"""
    f = TimeFunction(func::Function)

1-argument function of time that can be evaluated as `f(t)`.

Can be used in functions like `interpolate!` and `integrate` to dispatch
on the intended input arguments.

# Example:
```julia-repl
julia> f = TimeFunction(t -> t > 10.0 ? 1 : 0);

julia> f(1.0)
0

julia> f(20.0)
1
```
"""
struct TimeFunction{F}
    f::F
end

(f::TimeFunction{F})(t) where {F} = f.f(t)
(f::TimeFunction{F})(t) where {F <: Number} = f.f
(f::TimeFunction{F})(t) where {F <: AbstractArray} = f.f
(f::TimeFunction{F})(t) where {F <: RefValue} = f.f[]

"""
    f = SpaceTimeFunction(func::Function)

2-argument function of space and time that can be evaluated as `f(x, t)`.

Can be used in functions like `interpolate!` and `integrate` to dispatch
on the intended input arguments.

# Example:
```julia-repl
julia> f = SpaceTimeFunction((x, t) -> x[1] + t)
SpaceTimeFunction{var"#7#8"}(var"#7#8"())

julia> f([1, 2, 3], 1.123)
2.123
```
"""
struct SpaceTimeFunction{F}
    f::F
end

(f::SpaceTimeFunction{F})(x, t) where {F} = f.f(x, t)
(f::SpaceTimeFunction{F})(x, t) where {F <: Number} = f.f
(f::SpaceTimeFunction{F})(x, t) where {F <: AbstractArray} = f.f
(f::SpaceTimeFunction{F})(x, t) where {F <: RefValue} = f.f[]

####################################
# Field Interpolation
####################################

# Type parameter `interpolant::T where T` to force specialization on argument, see:
# https://docs.julialang.org/en/v1/manual/performance-tips/#Be-aware-of-when-Julia-avoids-specializing
# Probably only relevant in top scope usage of `interpolate!` (e.g. benchmarking of `interpolate!`).

"""
    interpolate!(field, interpolant)
    interpolate!(field, interpolant, tags)

Interpolates the the `field` using the given `interpolant`. `tags` can be passen optionally
to specify onto which parts of the field the interpolation should be applied.

`interpolant` can be a `Number`, `AbstractArray`, `RefValue` or `Function`
or a user-defined type that implements the method `eval_interpolant`.

If `interpolant` is a `Function`, its signature must be

    (field, proxy, q) -> value

where `q` is the index of the evaluation point.

Interpolant can also be a `SpaceFunction`, `TimeFunction` or `SpaceTimeFunction`.
"""
function interpolate!(field::AbstractDomainField, interpolant::T, tags=field.tags) where {T}
    fill!(field.assigned, false)
    for index in eachelementindex(field.mesh, tags)
        setproxy_dof_evaluation!(field, index)
        proxy = field.proxies[threadid()]
        for (q, dof) in enumerate(proxy.dofs)
            field.assigned[dof] == true && continue
            field.assigned[dof] = true
            field.u[dof] = eval_interpolant(field, interpolant, q)
        end
    end
    return field
end

function interpolate!(field::AbstractBoundaryField, interpolant::T, tags=field.tags) where {T}
    fill!(field.assigned, false)
    for index in eachfacetindex(field.mesh, tags)
        setproxy_dof_evaluation!(field, index)
        proxy = field.proxies[threadid()]
        for (q, dof) in enumerate(proxy.dofs)
            field.assigned[dof] == true && continue
            field.assigned[dof] = true
            field.u[dof] = eval_interpolant(field, interpolant, q)
        end
    end
    return field
end

"""
    interpolate_boundary!(field, interpolant, tags)

Interpolates the boundaries of the `field` specified by `tags`
using the given `interpolant`.

`interpolant` can be a `Number`, `AbstractArray`, `RefValue` or `Function`
or a user-defined type that implements the method `eval_interpolant`.

If `interpolant` is a `Function`, its signature must be

    (x, t) -> value
"""
function interpolate_boundary!(field::AbstractDomainField, interpolant::T, tags) where {T}
    fill!(field.assigned, false)
    for index in eachfacetindex(field.mesh, tags)
        setproxy_dof_evaluation!(field, index)
        proxy = field.proxies[threadid()]
        for (q, dof) in enumerate(view(proxy.dofs, proxy.localdofs))
            field.assigned[dof] == true && continue
            field.assigned[dof] = true
            field.u[dof] = eval_interpolant(field, interpolant, q)
        end
    end
    return field
end

function interpolate_boundary!(field::AbstractBoundaryField, args...; kwargs...)
    interpolate!(field, args...; kwargs...)
end

####################################
# Interpolant Evaluation
####################################

function eval_interpolant(field::AbstractField, interpolant, q::Integer)
    proxy = field.proxies[threadid()]
    eval_interpolant(field, proxy, interpolant, q)
end

function eval_interpolant(field::AbstractField, proxy::AbstractProxy, interpolant::Union{Number, AbstractArray}, q::Integer)
    return interpolant
end

function eval_interpolant(field::AbstractField, proxy::AbstractProxy, interpolant::RefValue, q::Integer)
    return interpolant[]
end

function eval_interpolant(field::AbstractField, proxy::AbstractProxy, interpolant::Function, q::Integer)
    return interpolant(field, proxy, q)
end

function eval_interpolant(field::AbstractField, proxy::AbstractProxy, interpolant::SpaceFunction, q::Integer)
    x = proxy.x[q]
    return interpolant(x)
end

function eval_interpolant(field::AbstractField, proxy::AbstractProxy, interpolant::TimeFunction, q::Integer)
    t = field.t[]
    return interpolant(t)
end

function eval_interpolant(field::AbstractField, proxy::AbstractProxy, interpolant::SpaceTimeFunction, q::Integer)
    x = proxy.x[q]
    t = field.t[]
    return interpolant(x, t)
end

####################################
# Field - Field Interpolation
####################################

defaultinterpolant(sourcefield, sourceproxy, q) = interpolate(sourcefield, q)
defaultinterpolant_grad(sourcefield, sourceproxy, q) = interpolate_grad(sourcefield, q)

"""
"""
function interpolate!(targetfield::AbstractDomainField, sourcefield::AbstractField, args...)
    interpolate!(defaultinterpolant, targetfield, sourcefield, args...)
end

function interpolate!(targetfield::AbstractBoundaryField, sourcefield::AbstractField, args...)
    interpolate!(defaultinterpolant, targetfield, sourcefield, args...)
end

"""
"""
function interpolate_grad!(targetfield::AbstractField, sourcefield::AbstractField, args...)
    interpolate_grad!(defaultinterpolant_grad, targetfield, sourcefield, args...)
end

function interpolate_grad!(interpolant::T, targetfield::AbstractField, sourcefield::AbstractField, args...) where {T}
    interpolate_nodalmean!(interpolant, targetfield, sourcefield, args...)
end

"""
"""
function interpolate!(interpolant::T, targetfield::AbstractDomainField, sourcefield::AbstractDomainField, tags=sourcefield.tags) where {T}
    fill!(targetfield.assigned, false)
    if first(targetfield.proxies) === first(sourcefield.proxies)
        for index in eachshapeindex(sourcefield, tags)
            setproxy_dof_evaluation!(targetfield, index)
            targetproxy = targetfield.proxies[threadid()]
            for (q, dof) in enumerate(targetproxy.dofs)
                targetfield.assigned[dof] == true && continue
                targetfield.assigned[dof] = true
                targetfield.u[dof] = eval_interpolant(sourcefield, interpolant, q)
            end
        end
    else
        for index in eachshapeindex(sourcefield, tags)
            setproxy_dof_evaluation!(targetfield, index)
            targetproxy = targetfield.proxies[threadid()]
            initproxy_point_evaluation!(sourcefield, index, first(targetproxy._evalpoints))
            for (q, dof) in enumerate(targetproxy.dofs)
                setproxy_point_evaluation!(sourcefield, index)
                targetfield.assigned[dof] == true && continue
                targetfield.assigned[dof] = true
                targetfield.u[dof] = eval_interpolant(sourcefield, interpolant, q)
            end
        end
    end
    return targetfield
end

function interpolate!(interpolant::T, targetfield::AbstractBoundaryField, sourcefield::AbstractBoundaryField, tags=sourcefield.tags) where {T}
    fill!(targetfield.assigned, false)
    if first(targetfield.proxies) === first(sourcefield.proxies)
        for index in eachshapeindex(sourcefield, tags)
            setproxy_dof_evaluation!(targetfield, index)
            targetproxy = targetfield.proxies[threadid()]
            for (q, dof) in enumerate(targetproxy.dofs)
                targetfield.assigned[dof] == true && continue
                targetfield.assigned[dof] = true
                targetfield.u[dof] = eval_interpolant(sourcefield, interpolant, q)
            end
        end
    else
        for index in eachshapeindex(sourcefield, tags)
            setproxy_dof_evaluation!(targetfield, index)
            targetproxy = targetfield.proxies[threadid()]
            initproxy_point_evaluation!(sourcefield, index, first(targetproxy._evalpoints))
            for (q, dof) in enumerate(targetproxy.dofs)
                setproxy_point_evaluation!(sourcefield, index)
                targetfield.assigned[dof] == true && continue
                targetfield.assigned[dof] = true
                targetfield.u[dof] = eval_interpolant(sourcefield, interpolant, q)
            end
        end
    end
    return targetfield
end

function interpolate!(interpolant::T, targetfield::AbstractDomainField, sourcefield::AbstractBoundaryField, tags=sourcefield.tags) where {T}
    fill!(targetfield.assigned, false)
    for index in eachshapeindex(sourcefield, tags)
        setproxy_dof_evaluation!(targetfield, index)
        targetproxy = targetfield.proxies[threadid()]
        Fs = typeof(sourcefield.mesh[index])
        for (q, localdof) in enumerate(targetproxy.localdofs)
            dof = targetproxy.dofs[localdof]
            ξ = dof_refcoordinates(Fs, targetfield.shapefunc, q)
            setproxy_point_evaluation!(sourcefield, index, ξ)
            targetfield.assigned[dof] == true && continue
            targetfield.assigned[dof] = true
            targetfield.u[dof] = eval_interpolant(sourcefield, interpolant, 1)
        end
    end
    return targetfield
end

function interpolate!(interpolant::T, targetfield::AbstractBoundaryField, sourcefield::AbstractDomainField, tags=targetfield.tags) where {T}
    fill!(targetfield.assigned, false)
    for index in eachshapeindex(targetfield, tags)
        setproxy_dof_evaluation!(targetfield, index)
        targetproxy = targetfield.proxies[threadid()]
        S = typeof(sourcefield.mesh[ElementIndex(index.elementid)])
        Fs = typeof(sourcefield.mesh[index])
        for (q, dof) in enumerate(targetproxy.dofs)
            _ξ = dof_refcoordinates(Fs, targetfield.shapefunc, q)
            ξ = facet2elementcoordinates(S, index.ilocal, _ξ)
            setproxy_point_evaluation!(sourcefield, ElementIndex(index.elementid), ξ)
            targetfield.assigned[dof] == true && continue
            targetfield.assigned[dof] = true
            targetfield.u[dof] = eval_interpolant(sourcefield, interpolant, 1)
        end
    end
    return targetfield
end

"""
"""
function interpolate_nodalmean!(targetfield::AbstractField, sourcefield::AbstractField, args...)
    interpolate_nodalmean!(defaultinterpolant, targetfield, sourcefield, args...)
end

function interpolate_nodalmean!(interpolant::T, targetfield::AbstractDomainField, sourcefield::AbstractDomainField, tags=sourcefield.tags) where {T}
    _assignnodalweights!(targetfield, eachshapeindex(sourcefield, tags))
    for i in eachindex(targetfield.u)
        if targetfield.assigned[i] > 0
            targetfield.u[i] = zero(eltype(targetfield.u))
        end
    end

    if first(targetfield.proxies) === first(sourcefield.proxies)
        for index in eachshapeindex(sourcefield, tags)
            setproxy_dof_evaluation!(targetfield, index)
            targetproxy = targetfield.proxies[threadid()]
            for (q, dof) in enumerate(targetproxy.dofs)
                targetfield.u[dof] += eval_interpolant(sourcefield, interpolant, q) / targetfield.assigned[dof]
            end
        end
    else
        for index in eachshapeindex(sourcefield, tags)
            setproxy_dof_evaluation!(targetfield, index)
            targetproxy = targetfield.proxies[threadid()]
            initproxy_point_evaluation!(sourcefield, index, first(targetproxy._evalpoints))
            for (q, dof) in enumerate(targetproxy.dofs)
                setproxy_point_evaluation!(sourcefield, index)
                targetfield.u[dof] += eval_interpolant(sourcefield, interpolant, q) / targetfield.assigned[dof]
            end
        end
    end
    return targetfield
end

function interpolate_nodalmean!(interpolant::T, targetfield::AbstractBoundaryField, sourcefield::AbstractBoundaryField, tags=sourcefield.tags) where {T}
    _assignnodalweights!(targetfield, eachshapeindex(sourcefield, tags))
    for i in eachindex(targetfield.u)
        if targetfield.assigned[i] > 0
            targetfield.u[i] = zero(eltype(targetfield.u))
        end
    end

    if first(targetfield.proxies) === first(sourcefield.proxies)
        for index in eachshapeindex(sourcefield, tags)
            setproxy_dof_evaluation!(targetfield, index)
            targetproxy = targetfield.proxies[threadid()]
            for (q, dof) in enumerate(targetproxy.dofs)
                targetfield.u[dof] += eval_interpolant(sourcefield, interpolant, q) / targetfield.assigned[dof]
            end
        end
    else
        for index in eachshapeindex(sourcefield, tags)
            setproxy_dof_evaluation!(targetfield, index)
            targetproxy = targetfield.proxies[threadid()]
            initproxy_point_evaluation!(sourcefield, index, first(targetproxy._evalpoints))
            for (q, dof) in enumerate(targetproxy.dofs)
                setproxy_point_evaluation!(sourcefield, index)
                targetfield.u[dof] += eval_interpolant(sourcefield, interpolant, q) / targetfield.assigned[dof]
            end
        end
    end
    return targetfield
end

function interpolate_nodalmean!(interpolant::T, targetfield::AbstractDomainField, sourcefield::AbstractBoundaryField, tags=sourcefield.tags) where {T}
    _assignnodalweights_facets!(targetfield, eachshapeindex(sourcefield, tags))
    for i in eachindex(targetfield.u)
        if targetfield.assigned[i] > 0
            targetfield.u[i] = zero(eltype(targetfield.u))
        end
    end

    for index in eachshapeindex(sourcefield, tags)
        setproxy_dof_evaluation!(targetfield, index)
        targetproxy = targetfield.proxies[threadid()]
        Fs = typeof(sourcefield.mesh[index])
        for (q, localdof) in enumerate(targetproxy.localdofs)
            dof = targetproxy.dofs[localdof]
            ξ = dof_refcoordinates(Fs, targetfield.shapefunc, q)
            setproxy_point_evaluation!(sourcefield, index, ξ)
            targetfield.u[dof] += eval_interpolant(sourcefield, interpolant, 1) / targetfield.assigned[dof]
        end
    end
    return targetfield
end

function interpolate_nodalmean!(interpolant::T, targetfield::AbstractBoundaryField, sourcefield::AbstractDomainField, tags=targetfield.tags) where {T}
    _assignnodalweights!(targetfield, eachshapeindex(targetfield, tags))
    for i in eachindex(targetfield.u)
        if targetfield.assigned[i] > 0
            targetfield.u[i] = zero(eltype(targetfield.u))
        end
    end

    for index in eachshapeindex(targetfield, tags)
        setproxy_dof_evaluation!(targetfield, index)
        targetproxy = targetfield.proxies[threadid()]
        S = typeof(sourcefield.mesh[ElementIndex(index.elementid)])
        Fs = typeof(sourcefield.mesh[index])
        for (q, dof) in enumerate(targetproxy.dofs)
            _ξ = dof_refcoordinates(Fs, targetfield.shapefunc, q)
            ξ = facet2elementcoordinates(S, index.ilocal, _ξ)
            setproxy_point_evaluation!(sourcefield, ElementIndex(index.elementid), ξ)
            targetfield.u[dof] += eval_interpolant(sourcefield, interpolant, 1) / targetfield.assigned[dof]
        end
    end
    return targetfield
end

function _assignnodalweights!(field, indices)
    fill!(field.assigned, 0)
    shapes = (field.mesh[index] for index in indices)
    for shape in shapes
        dofs = field.dofs[shape.id]
        for dof in dofs
            field.assigned[dof] += 1
        end
    end
    return field.assigned
end

function _assignnodalweights_facets!(field, facetindices)
    fill!(field.assigned, 0)
    for index in facetindices
        dofs = field.dofs[index.elementid]
        for dof in dofs
            field.assigned[dof] += 1
        end
    end
    return field.assigned
end
