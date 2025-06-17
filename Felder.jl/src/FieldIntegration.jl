export integrate
export integrate_boundary
export integrate_edge # TODO

# TODO: Checks that the tags actually exists in the field

defaultintegrand(field, proxy, q) = interpolate(field, q)

# Type parameter `interpolant::T where T` to force specialization on argument, see:
# https://docs.julialang.org/en/v1/manual/performance-tips/#Be-aware-of-when-Julia-avoids-specializing
# Probably only relevant in top scope usage of `interpolate!` (e.g. benchmarking of `interpolate!`).

"""
    integrate(field)
    integrate(field, tags)
    integrate(interpolant, field, tags)
    integrate(interpolant, field, tags, quadrature)

Integrates the `integrand` over the domain of `field`.

`integrand` can be a `Number`, `AbstractArray`, `RefValue` or `Function`
or a user-defined type that implements the method `eval_interpolant`.

The default integrand is the field variable `field.u`.

If `integrand` is a `Function`, its signature must be

    (field, proxy, t) -> value.

where `q` is the index of the evaluation point.

Set the `integrand` to 1.0 to calculate the volume, area, or length of a
domain, boundary or edge, respectively.

Supports do-block syntax, e.g.

    integrate(field) do field, proxy, q
        # do something, e.g. return proxy.x[q] * interpolate(field, q)
    end
"""
function integrate(field::AbstractField, args...)
    integrate(defaultintegrand, field, args...)
end

function integrate(integrand::T, field::AbstractDomainField, tags=field.tags,
        quadrature=defaultquadrature(field)) where {T}

    val = zero(typeof_integrand(field, integrand))
    for index in eachelementindex(field.mesh, tags)
        setproxy_integration!(field, index)
        proxy = field.proxies[threadid()]
        for q in eachindex(proxy.dΩ)
            val += eval_integrand(field, integrand, q) * proxy.dΩ[q]
        end
    end
    return val
end

function integrate(integrand::T, field::AbstractBoundaryField, tags=field.tags,
        quadrature=defaultquadrature(field)) where {T}

    val = zero(typeof_integrand(field, integrand))
    for index in eachfacetindex(field.mesh, tags)
        setproxy_integration!(field, index)
        proxy = field.proxies[threadid()]
        for q in eachindex(proxy.dΩ)
            val += eval_integrand(field, integrand, q) * proxy.dΩ[q]
        end
    end
    return val
end


"""
    integrate_boundary(field, tags)
    integrate_boundary(integrand, field, tags)
    integrate_boundary(integrand, field, tags, quadrature)

Integrates the `integrand` over the boundaries of `field`.

`integrand` can be a `Number`, `AbstractArray`, `RefValue` or `Function`
or a user-defined type that implements the method `eval_interpolant`.

The default integrand is the field variable `field.u`.

If `integrand` is a `Function`, its signature must be

    (field, proxy, t) -> value.

Supports do-block syntax, e.g.

    integrate_boundary(field) do field, proxy, q
        # do something, e.g. return proxy.normals[q] * interpolate(field, q)
    end
"""
function integrate_boundary(field::AbstractField, args...)
    integrate_boundary((field, proxy, q) -> interpolate(field, q), field, args...)
end

function integrate_boundary(integrand::T, field::AbstractDomainField, tags,
        quadrature=defaultquadrature(field)) where {T}

    val = zero(typeof_integrand(field, integrand))
    for index in eachfacetindex(field.mesh, tags)
        setproxy_integration!(field, index)
        proxy = field.proxies[threadid()]
        for q in eachindex(proxy.dΩ)
            val += eval_integrand(field, integrand, q) * proxy.dΩ[q]
        end
    end
    return val
end

#-------------------------

function typeof_integrand(field::AbstractDomainField, integrand)
    setproxy_integration!(field, first(eachelementindex(field.mesh, field.tags)))
    mockval = eval_integrand(field, integrand, 1)
    return typeof(mockval)
end

function typeof_integrand(field::AbstractBoundaryField, integrand)
    setproxy_integration!(field, first(eachfacetindex(field.mesh, field.tags)))
    mockval = eval_integrand(field, integrand, 1)
    return typeof(mockval)
end

#-------------------------

function eval_integrand(field::AbstractField, integrand, q::Integer)
    proxy = field.proxies[threadid()]
    return eval_integrand(field, proxy, integrand, q)
end

function eval_integrand(field::AbstractField, proxy::AbstractProxy, integrand::Union{Number, AbstractArray}, q::Integer)
    return integrand
end

function eval_integrand(field::AbstractField, proxy::AbstractProxy, integrand::RefValue, q::Integer)
    return integrand[]
end

function eval_integrand(field::AbstractField, proxy::AbstractProxy, integrand::Function, q::Integer)
    return integrand(field, proxy, q)
end
