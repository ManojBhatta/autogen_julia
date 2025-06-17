export AbstractProxy
export Proxy
export Proxy1D
export Proxy2D
export Proxy3D
export ScalarProxy
export ScalarProxy1D
export ScalarProxy2D
export ScalarProxy3D
export setdofs!
export set_element_integration!
export set_facet_integration!
export set_element_dof_evaluation!
export set_facet_dof_evaluation!
export set_point_evaluation!
export init_point_evaluation!
export interpolate
export interpolate_cache
export interpolate_grad
export interpolate_comp_grad
export nrefdims

#############################
# Proxy
#############################

abstract type AbstractProxy{Dim, RefDim, Sf} end

Base.ndims(::AbstractProxy{Dim}) where {Dim} = Dim
Base.ndims(::Type{<:AbstractProxy{Dim}}) where {Dim} = Dim
nrefdims(::Type{<:AbstractProxy{Dim, RefDim}}) where {Dim, RefDim} = RefDim

sftype(::AbstractProxy{Dim, RefDim, Sf}) where {Dim, RefDim, Sf} = Sf
sftype(::Type{<:AbstractProxy{Dim, RefDim, Sf}}) where {Dim, RefDim, Sf} = Sf

ndofs(proxy::AbstractProxy) = length(proxy.dofs)

function setdofs!(proxy::AbstractProxy, dofs)
    # sortperm allocates for Julia < 1.9! https://github.com/JuliaLang/julia/pull/47966
    # allocates when alg= is not given in Julia 1.9. TODO: Check performance of insertion sort
    proxy.dofs .= dofs
    sortperm!(proxy.dofperm, dofs, alg=InsertionSort)
    return proxy
end

function setdofs_noperm!(proxy::AbstractProxy, dofs)
    proxy.dofs .= dofs
    return proxy
end

abstract type AbstractScalarProxy{Dim, RefDim, Sf} <: AbstractProxy{Dim, RefDim, Sf} end
abstract type AbstractVectorProxy{Dim, RefDim, Sf} <: AbstractProxy{Dim, RefDim, Sf} end

"""
"""
struct ScalarProxy{Dim, RefDim, Sf} <: AbstractScalarProxy{Dim, RefDim, Sf}
    shapefunc::Sf
    _shapetype::RefValue{Int}
    _mode::RefValue{Int} # 0 for element integration; 1, 2, 3, ... for facets integration; -1, -2, ... for lagrange facet evaluation; 10 for lagrange element evaluation, 20 for multi point evaluation, 21 for single point evaluation
    _quadorder::RefValue{Int}
    _shapeid::RefValue{Int}

    _evalpoints::Vector{Vector{SVector{RefDim, Float64}}}
    _quadweights::Vector{Vector{Float64}}

    _N::Vector{Vector{Vector{Float64}}}
    _dNdξ::Vector{Vector{Vector{SVector{RefDim, Float64}}}}

    _M::Vector{Vector{Vector{Float64}}}
    _dMdξ::Vector{Vector{Vector{SVector{RefDim, Float64}}}}

    x::Vector{SVector{Dim, Float64}}
    N::Vector{Vector{Float64}}
    dNdx::Vector{Vector{SVector{Dim, Float64}}}
    dΩ::Vector{Float64}
    normals::Vector{SVector{Dim, Float64}}

    dofs::Vector{Int}
    dofperm::Vector{Int}

    _localdofs::Vector{Vector{Int}}
    localdofs::Vector{Int}
end

const ScalarProxy1D = ScalarProxy{1}
const ScalarProxy2D = ScalarProxy{2}
const ScalarProxy3D = ScalarProxy{3}

const Proxy = ScalarProxy
const Proxy1D = Proxy{1}
const Proxy2D = Proxy{2}
const Proxy3D = Proxy{3}

"""
"""
function ScalarProxy{Dim, RefDim, Sf}(shapefunc::Sf) where {Dim, RefDim, Sf <: AbstractScalarShapeFunctions}
    _shapetype = Ref(0)
    _mode = Ref(0)
    _quadorder = Ref(0)
    _shapeid = Ref(0)

    _evalpoints  = Vector{Vector{SVector{RefDim, Float64}}}()
    _quadweights = Vector{Vector{Float64}}()

    _N   = Vector{Vector{Vector{Float64}}}()
    _dNdξ = Vector{Vector{Vector{SVector{RefDim, Float64}}}}()

    _M   = Vector{Vector{Vector{Float64}}}()
    _dMdξ = Vector{Vector{Vector{SVector{RefDim, Float64}}}}()

    x    = Vector{SVector{Dim, Float64}}()
    N    = Vector{Vector{Float64}}()
    dNdx = Vector{Vector{SVector{Dim, Float64}}}()
    dΩ   = Vector{Float64}()
    normals = Vector{SVector{Dim, Float64}}()

    dofs = Vector{Int}()
    dofperm = Vector{Int}()

    _localdofs = Vector{Vector{Int}}()
    localdofs = Vector{Int}()

    ScalarProxy{Dim, RefDim, typeof(shapefunc)}(
        shapefunc,
        _shapetype,
        _mode,
        _quadorder,
        _shapeid,
        _evalpoints,
        _quadweights,
        _N,
        _dNdξ,
        _M,
        _dMdξ,
        x,
        N,
        dNdx,
        dΩ,
        normals,
        dofs,
        dofperm,
        _localdofs,
        localdofs,
    )
end

function ScalarProxy{Dim}(args...; kwargs...) where {Dim}
    ScalarProxy{Dim, Dim}(args...; kwargs...)
end

function ScalarProxy{Dim, RefDim}(shapefunc::Sf, args...; kwargs...) where {Dim, RefDim, Sf}
    ScalarProxy{Dim, RefDim, Sf}(shapefunc, args...; kwargs...)
end

function Base.show(io::IO, proxy::ScalarProxy{Dim, RefDim, Sf}) where {Dim, RefDim, Sf}
    # TODO
    println(io, "ScalarProxy{$Dim, $RefDim, $Sf}", " with")
    println(io, "  Proxy mode:                      ", proxy._mode[])
    println(io, "  Number of geometric nodes:       ", proxy._shapetype[])
    println(io, "  Quadrature order:                ", proxy._quadorder[])
    println(io, "  Shape ID:                        ", proxy._shapeid[])
    println(io, "  Number of degrees of freedom:    ", length(proxy.dofs))
    println(io, "  Degrees of freedom:              ", proxy.dofs)
    println(io, "  Number of eval. point sets:      ", length(proxy._evalpoints))
    for i in 1:length(proxy._evalpoints)
       println(io,  "    Set $i: Number of eval. points: ", length(proxy._evalpoints[i]))
    end
end

function resize_proxy!(proxy::AbstractScalarProxy, n_sets, n_evalpoints, n_dofs, n_dofs_geo, n_localdofs=0)
    nested_resize_nan!(proxy._evalpoints, n_sets, n_evalpoints)
    nested_resize_nan!(proxy._quadweights, n_sets, n_evalpoints)

    nested_resize_nan!(proxy._N, n_sets, n_evalpoints, n_dofs)
    nested_resize_nan!(proxy._dNdξ, n_sets, n_evalpoints, n_dofs)

    nested_resize_nan!(proxy._M, n_sets, n_evalpoints, n_dofs_geo)
    nested_resize_nan!(proxy._dMdξ, n_sets, n_evalpoints, n_dofs_geo)

    resize_nan!(proxy.x, n_evalpoints)
    nested_resize_nan!(proxy.N, n_evalpoints, n_dofs)
    nested_resize_nan!(proxy.dNdx, n_evalpoints, n_dofs)
    resize_nan!(proxy.dΩ, n_evalpoints)
    resize_nan!(proxy.normals, n_evalpoints)

    resize_zero!(proxy.dofs, n_dofs)
    resize_zero!(proxy.dofperm, n_dofs)

    nested_resize_zero!(proxy._localdofs, n_sets, n_localdofs)
    resize_zero!(proxy.localdofs, n_localdofs)

    return proxy
end

function _to_globalcoordinates(shapecoordinates, shapefunctions_M)
    @boundscheck @assert length(shapecoordinates) == length(shapefunctions_M)  # TODO boundscheck
    x = zero(first(shapecoordinates))
     for (xi, Ni) in zip(shapecoordinates, shapefunctions_M)
        @inbounds x += xi * Ni
    end
    return x
end

function _parameterization_jacobian(shapecoordinates, shapefunctions_∇M)
    @boundscheck @assert length(shapecoordinates) == length(shapefunctions_∇M) # TODO boundscheck
    N = length(first(shapecoordinates))
    M = length(first(shapefunctions_∇M))
    J = zero(SMatrix{N, M, Float64})
    for (xi, dNdξi) in zip(shapecoordinates, shapefunctions_∇M)
        @inbounds J += xi .* dNdξi'
    end
    return J
end

############################
# Element Integration
############################

function set_element_integration!(proxy::AbstractProxy, ::Type{S}, shapecoordinates,
        quadrature=defaultquadrature(proxy.shapefunc, S)) where {S<:AbstractShape}

    _init_element_integration!(proxy, S, quadrature)
    _set_element_integration!(proxy, S, shapecoordinates, quadrature)
end

function _init_element_integration!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S},
        quadrature::AbstractQuadrature) where {Dim, RefDim, S<:AbstractShape{RefDim}}

    if proxy._shapetype[] == nnodes(S) && proxy._mode[] == 0 && proxy._quadorder[] == getorder(quadrature)
        return proxy
    end

    sf = proxy.shapefunc
    proxy._shapetype[] = nnodes(S)
    proxy._mode[] = 999
    proxy._quadorder[] = getorder(quadrature)

    resize_proxy!(proxy, 1, nquadpoints(S, quadrature), ndofs(S, sf), ndofs(S))

    for (q, ξ) in enumerate(eachquadpoint(S, quadrature)) # TODO: inbounds
        proxy._evalpoints[1][q] = ξ
        proxy._quadweights[1][q] = quadweight(S, quadrature, q)

        for i in 1:ndofs(S, sf)
            proxy._N[1][q][i] = shapefunc_N(S, sf, i, ξ)
            proxy._dNdξ[1][q][i] = shapefunc_∇N(S, sf, i, ξ)
        end

        for i in 1:ndofs(S)
            proxy._M[1][q][i] = shapefunc_N(S, sftype(S)(), i, ξ)
            proxy._dMdξ[1][q][i] = shapefunc_∇N(S, sftype(S)(), i, ξ)
        end
    end

    proxy.N .= proxy._N[1]

    return proxy
end

function _set_element_integration!(proxy::AbstractScalarProxy{Dim, Dim}, ::Type{S}, shapecoordinates,
        quadrature=defaultquadrature(proxy.shapefunc, S)) where {Dim, S<:AbstractShape{Dim}}

    @assert nnodes(S) == length(shapecoordinates)

    _evalpoints, _quadweights = proxy._evalpoints[1], proxy._quadweights[1]
    _dNdξ, _M, _dMdξ = proxy._dNdξ[1], proxy._M[1], proxy._dMdξ[1]

    for q in eachindex(_evalpoints) # TODO: inbounds
        proxy.x[q] = _to_globalcoordinates(shapecoordinates, _M[q])

        J = _parameterization_jacobian(shapecoordinates, _dMdξ[q])
        invJ = inv(J)'

        for i in eachindex(proxy.dNdx[q])
            proxy.dNdx[q][i] = invJ * _dNdξ[q][i]
        end

        proxy.dΩ[q] = _quadweights[q] * abs(det(J))
    end

    proxy._mode[] = 0

    return proxy
end

function _set_element_integration!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S}, shapecoordinates,
        quadrature=defaultquadrature(proxy.shapefunc, S)) where {Dim, RefDim, S<:AbstractShape{RefDim}}

    @assert nnodes(S) == length(shapecoordinates)

    _evalpoints, _quadweights = proxy._evalpoints[1], proxy._quadweights[1]
    _M, _dMdξ = proxy._M[1], proxy._dMdξ[1]

    for q in eachindex(_evalpoints) # TODO: inbounds
        proxy.x[q] = _to_globalcoordinates(shapecoordinates, _M[q])

        J = _parameterization_jacobian(shapecoordinates, _dMdξ[q])
        drds_x_drdt = parameterization_normal(S, J)
        norm_drds_x_drdt = norm(drds_x_drdt)

        proxy.normals[q] = drds_x_drdt / norm_drds_x_drdt
        proxy.dΩ[q] = _quadweights[q] * norm_drds_x_drdt
    end

    proxy._mode[] = 0

    return proxy
end

############################
# Facet Integration
############################

function set_facet_integration!(proxy::AbstractProxy, ::Type{S}, shapecoordinates, ifacet,
        quadrature=defaultquadrature(proxy.shapefunc, S)) where {S<:AbstractShape}

    _init_facet_integration!(proxy, S, quadrature)
    _set_facet_integration!(proxy, S, shapecoordinates, ifacet, quadrature)
end

function _init_facet_integration!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S},
        quadrature::AbstractQuadrature) where {Dim, RefDim, S<:AbstractShape{RefDim}}

    if proxy._shapetype[] == nnodes(S) && (0 < proxy._mode[] < 10) && proxy._quadorder[] == getorder(quadrature)
        return proxy
    end

    sf = proxy.shapefunc
    proxy._shapetype[] = nnodes(S)
    proxy._mode[] = 999
    proxy._quadorder[] = getorder(quadrature)

    # All shapes have the same shape type on all facets, except Wedge and Pyramid
    n_evalpoints = maximum(nquadpoints(localfacettype(S, i), quadrature) for i in 1:nfacets(S))
    n_localdofs = maximum(ndofs(localfacettype(S, i), sf) for i in 1:nfacets(S))
    resize_proxy!(proxy, nfacets(S), n_evalpoints, ndofs(S, sf), ndofs(S), n_localdofs)

    for f in 1:nfacets(S)
        Fs = localfacettype(S, f)
        resize!(proxy._evalpoints[f], nquadpoints(Fs, quadrature))
        resize!(proxy._quadweights[f], nquadpoints(Fs, quadrature))
        resize!(proxy._N[f], nquadpoints(Fs, quadrature))
        resize!(proxy._dNdξ[f], nquadpoints(Fs, quadrature))
        resize!(proxy._M[f], nquadpoints(Fs, quadrature))
        resize!(proxy._dMdξ[f], nquadpoints(Fs, quadrature))
        for (q, _ξ) in enumerate(eachquadpoint(Fs, quadrature)) # TODO: inbounds
            ξ = facet2elementcoordinates(S, f, _ξ)
            proxy._evalpoints[f][q] = ξ
            proxy._quadweights[f][q] = quadweight(Fs, quadrature, q)

            for i in 1:ndofs(S, sf)
                proxy._N[f][q][i] = shapefunc_N(S, sf, i, ξ)
                proxy._dNdξ[f][q][i] = shapefunc_∇N(S, sf, i, ξ)
            end

            for i in 1:ndofs(S)
                proxy._M[f][q][i] = shapefunc_N(S, sftype(S)(), i, ξ)
                proxy._dMdξ[f][q][i] = shapefunc_∇N(S, sftype(S)(), i, ξ)
            end
        end
        resize!(proxy._localdofs[f], ndofs(Fs, sf))
        proxy._localdofs[f] .= facetdofs(S, sf, f)
    end

    return proxy
end

function _set_facet_integration!(proxy::AbstractScalarProxy{Dim, Dim}, ::Type{S}, shapecoordinates, ifacet,
        quadrature=defaultquadrature(proxy.shapefunc, S)) where {Dim, S<:AbstractShape{Dim}}

    @assert nnodes(S) == length(shapecoordinates)
    @assert 1 <= ifacet <= nfacets(S)

    _evalpoints, _quadweights = proxy._evalpoints[ifacet], proxy._quadweights[ifacet]
    _N, _dNdξ, _M, _dMdξ = proxy._N[ifacet], proxy._dNdξ[ifacet], proxy._M[ifacet], proxy._dMdξ[ifacet]

    # Vectors of vectors like proxy.x, .N, .dNdx and .normals are not resized
    # to avoid allocating new vector references (when growing).
    # Therefore, for Wedges where the top and bottom facets have different number
    # of integration points, only proxy.dΩ has the correct size.
    resize!(proxy.dΩ, length(_evalpoints))

    for q in eachindex(_evalpoints) # TODO: inbounds
        proxy.x[q] = _to_globalcoordinates(shapecoordinates, _M[q])

        J = _parameterization_jacobian(shapecoordinates, _dMdξ[q])
        invJ = inv(J)'

        drds_x_drdt = parameterization_normal(S, J, ifacet)
        norm_drds_x_drdt = norm(drds_x_drdt)

        for i in eachindex(proxy.dNdx[q])
            proxy.dNdx[q][i] = invJ * _dNdξ[q][i]
        end

        proxy.normals[q] = drds_x_drdt / norm_drds_x_drdt
        proxy.dΩ[q] = _quadweights[q] * norm_drds_x_drdt
        proxy.N[q] = _N[q]
    end

    resize!(proxy.localdofs, length(proxy._localdofs[ifacet]))
    proxy.localdofs .= proxy._localdofs[ifacet]

    proxy._mode[] = ifacet

    return proxy
end

####################################
# DOF Point Element Evaluation
####################################
# normals are not set because the Lagrange evaluation points are not located on facets in general.
# Some are, but some points are located on the nodes, edge (3D) or the interior of the element (higher order)
# where no normal vector is defined.
# _quadweights, dΩ are not set because no integration is intended.

function set_element_dof_evaluation!(proxy::AbstractProxy, ::Type{S}, shapecoordinates) where {S<:AbstractShape}
    _init_element_dof_evaluation!(proxy, S)
    _set_element_dof_evaluation!(proxy, S, shapecoordinates)
end

function _init_element_dof_evaluation!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S}) where {Dim, RefDim, S<:AbstractShape{RefDim}}
    if proxy._shapetype[] == nnodes(S) && proxy._mode[] == 10
        return proxy
    end

    sf = proxy.shapefunc
    proxy._shapetype[] = nnodes(S)
    proxy._mode[] = 999
    proxy._quadorder[] = 0

    resize_proxy!(proxy, 1, ndofs(S, sf), ndofs(S, sf), ndofs(S))

    for (q, ξ) in enumerate(each_dof_refcoordinates(S, sf)) # TODO: inbounds
        proxy._evalpoints[1][q] = ξ

        for i in 1:ndofs(S, sf)
            proxy._N[1][q][i] = shapefunc_N(S, sf, i, ξ)
            proxy._dNdξ[1][q][i] = shapefunc_∇N(S, sf, i, ξ)
        end

        for i in 1:ndofs(S)
            proxy._M[1][q][i] = shapefunc_N(S, sftype(S)(), i, ξ)
            proxy._dMdξ[1][q][i] = shapefunc_∇N(S, sftype(S)(), i, ξ)
        end
    end

    proxy.N .= proxy._N[1]

    return proxy
end

function _set_element_dof_evaluation!(proxy::AbstractScalarProxy{Dim, Dim}, ::Type{S}, shapecoordinates) where {Dim, S<:AbstractShape{Dim}}
    @assert nnodes(S) == length(shapecoordinates)

    _evalpoints = proxy._evalpoints[1]
    _dNdξ, _M, _dMdξ = proxy._dNdξ[1], proxy._M[1], proxy._dMdξ[1]

    for q in eachindex(_evalpoints) # TODO: inbounds
        proxy.x[q] = _to_globalcoordinates(shapecoordinates, _M[q])

        J = _parameterization_jacobian(shapecoordinates, _dMdξ[q])
        invJ = inv(J)'

        for i in eachindex(proxy.dNdx[q])
            proxy.dNdx[q][i] = invJ * _dNdξ[q][i]
        end
    end

    proxy._mode[] = 10

    return proxy
end

function _set_element_dof_evaluation!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S}, shapecoordinates) where {Dim, RefDim, S<:AbstractShape{RefDim}}
    @assert nnodes(S) == length(shapecoordinates)

    _evalpoints = proxy._evalpoints[1]
    _M, _dMdξ = proxy._M[1], proxy._dMdξ[1]

    for q in eachindex(_evalpoints) # TODO: inbounds
        proxy.x[q] = _to_globalcoordinates(shapecoordinates, _M[q])

        J = _parameterization_jacobian(shapecoordinates, _dMdξ[q])
        drds_x_drdt = parameterization_normal(S, J)
        norm_drds_x_drdt = norm(drds_x_drdt)

        proxy.normals[q] = drds_x_drdt / norm_drds_x_drdt
    end

    proxy._mode[] = 10

    return proxy
end

####################################
# Dof Point Facet Evaluation
####################################
# _quadweights, dΩ are not set because no integration is intended.

function set_facet_dof_evaluation!(proxy::AbstractProxy, ::Type{S}, shapecoordinates, ifacet) where {S<:AbstractShape}
    _init_facet_dof_evaluation!(proxy, S)
    _set_facet_dof_evaluation!(proxy, S, shapecoordinates, ifacet)
end

function _init_facet_dof_evaluation!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S}) where {Dim, RefDim, S<:AbstractShape{RefDim}}
    if proxy._shapetype[] == nnodes(S) && (-10 < proxy._mode[] < 0)
        return proxy
    end

    sf = proxy.shapefunc
    proxy._shapetype[] = nnodes(S)
    proxy._quadorder[] = 0

    # All shapes have the same shape type on all facets, except Wedge and Pyramid
    n_evalpoints = maximum(ndofs(localfacettype(S, i), sf) for i in 1:nfacets(S))
    resize_proxy!(proxy, nfacets(S), n_evalpoints, ndofs(S, sf), ndofs(S), n_evalpoints)

    for f in 1:nfacets(S)
        Fs = localfacettype(S, f)
        resize!(proxy._evalpoints[f], ndofs(Fs, sf))
        for (q, _ξ) in enumerate(each_dof_refcoordinates(Fs, sf)) # TODO: inbounds
            ξ = facet2elementcoordinates(S, f, _ξ)
            proxy._evalpoints[f][q] = ξ

            for i in 1:ndofs(S, sf)
                proxy._N[f][q][i] = shapefunc_N(S, sf, i, ξ)
                proxy._dNdξ[f][q][i] = shapefunc_∇N(S, sf, i, ξ)
            end

            for i in 1:ndofs(S)
                proxy._M[f][q][i] = shapefunc_N(S, sftype(S)(), i, ξ)
                proxy._dMdξ[f][q][i] = shapefunc_∇N(S, sftype(S)(), i, ξ)
            end
        end
        resize!(proxy._localdofs[f], ndofs(Fs, sf))
        proxy._localdofs[f] .= facetdofs(S, sf, f)
    end

    return proxy
end

function _set_facet_dof_evaluation!(proxy::AbstractScalarProxy{Dim, Dim}, ::Type{S}, shapecoordinates, ifacet,
        shapeid=0) where {Dim, S<:AbstractShape{Dim}}

    @assert nnodes(S) == length(shapecoordinates)
    @assert 1 <= ifacet <= nfacets(S)

    _evalpoints = proxy._evalpoints[ifacet]
    _dNdξ, _M, _dMdξ = proxy._dNdξ[ifacet], proxy._M[ifacet], proxy._dMdξ[ifacet]

    for q in eachindex(_evalpoints) # TODO: inbounds
        proxy.x[q] = _to_globalcoordinates(shapecoordinates, _M[q])

        J = _parameterization_jacobian(shapecoordinates, _dMdξ[q])
        invJ = inv(J)'

        drds_x_drdt = parameterization_normal(S, J, ifacet)
        norm_drds_x_drdt = norm(drds_x_drdt)

        for i in eachindex(proxy.dNdx[q])
            proxy.dNdx[q][i] = invJ * _dNdξ[q][i]
        end

        proxy.normals[q] = drds_x_drdt / norm_drds_x_drdt
    end

    proxy.N .= proxy._N[ifacet]

    resize!(proxy.localdofs, length(proxy._localdofs[ifacet]))
    proxy.localdofs .= proxy._localdofs[ifacet]

    proxy._mode[] = -ifacet

    return proxy
end

############################
# Multi Point Evaluation
############################
# The points ξs must be initialized with `init_point_evaluation!` once ]
# before calling `set_point_evaluation!` iteratively for all shapes.
# The points ξ cannot be changed when calling `set_point_evaluation!`,
# and the shape type must not change (otherwise error).
# _quadweights, dΩ are not set because no integration is intended.

function init_point_evaluation!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S}, ξs::AbstractVector{<:AbstractVector}) where {Dim, RefDim, S<:AbstractShape{RefDim}}
    if proxy._shapetype[] == nnodes(S) && proxy._mode[] == 20 && proxy._evalpoints[1] == ξs
        return proxy
    end

    @assert RefDim == length(first(ξs))

    sf = proxy.shapefunc
    proxy._shapetype[] = nnodes(S)
    proxy._mode[] = 20
    proxy._quadorder[] = 0
    proxy._shapeid[] = 0

    resize_proxy!(proxy, 1, length(ξs), ndofs(S, sf), ndofs(S))

    for (q, ξ) in enumerate(ξs) # TODO: inbounds
        proxy._evalpoints[1][q] = ξ

        for i in 1:ndofs(S, sf)
            proxy._N[1][q][i] = shapefunc_N(S, sf, i, ξ)
            proxy._dNdξ[1][q][i] = shapefunc_∇N(S, sf, i, ξ)
        end

        for i in 1:ndofs(S)
            proxy._M[1][q][i] = shapefunc_N(S, sftype(S)(), i, ξ)
            proxy._dMdξ[1][q][i] = shapefunc_∇N(S, sftype(S)(), i, ξ)
        end
    end

    proxy.N .= proxy._N[1]

    return proxy
end

# TODO: Check allocations and @btime
# Normal not evaluated because we don't know where the point is located (node, edge, interior)
function set_point_evaluation!(proxy::AbstractScalarProxy{Dim, Dim}, ::Type{S}, shapecoordinates) where {Dim, S<:AbstractShape{Dim}}
    @assert nnodes(S) == length(shapecoordinates)
    @assert proxy._mode[] == 20 "`proxy` must be initialized with `init_point_evaluation!`
        for multi point evaluation before calling `set_point_evaluation!`"
    @assert proxy._shapetype[] == nnodes(S) "Shape type cannot change when calling
        `set_point_evaluation!` for multi point evaluation"

    _evalpoints = proxy._evalpoints[1]
    _dNdξ, _M, _dMdξ = proxy._dNdξ[1], proxy._M[1], proxy._dMdξ[1]

    for q in eachindex(_evalpoints) # TODO: inbounds
        proxy.x[q] = _to_globalcoordinates(shapecoordinates, _M[q])

        J = _parameterization_jacobian(shapecoordinates, _dMdξ[q])
        invJ = inv(J)'

        for i in eachindex(proxy.dNdx[q])
            proxy.dNdx[q][i] = invJ * _dNdξ[q][i]
        end
    end

    return proxy
end

# Normal is evaluated because we are on a facet in higher dimensional space
function set_point_evaluation!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S}, shapecoordinates) where {Dim, RefDim, S<:AbstractShape{RefDim}}
    @assert nnodes(S) == length(shapecoordinates)
    @assert proxy._mode[] == 20 "`proxy` must be initialized with `init_point_evaluation!`
        for multi point evaluation before calling `set_point_evaluation!`"
    @assert proxy._shapetype[] == nnodes(S) "Shape type cannot change when calling
        `set_point_evaluation!` for multi point evaluation"

    _evalpoints = proxy._evalpoints[1]
    _M, _dMdξ = proxy._M[1], proxy._dMdξ[1]

    for q in eachindex(_evalpoints) # TODO: inbounds
        proxy.x[q] = _to_globalcoordinates(shapecoordinates, _M[q])

        J = _parameterization_jacobian(shapecoordinates, _dMdξ[q])
        drds_x_drdt = parameterization_normal(S, J)
        norm_drds_x_drdt = norm(drds_x_drdt)

        proxy.normals[q] = drds_x_drdt / norm_drds_x_drdt
    end

    return proxy
end

############################
# Single Point Evaluation
############################
# The point ξ and the shape type can change with each call to set_point_evaluation!.
# _quadweights, dΩ are not set because no integration is intended.
# Shape functions are evaluated for each call set_point_evaluation! with respect to ξ

function set_point_evaluation!(proxy::AbstractProxy, ::Type{S}, shapecoordinates, ξ::AbstractVector{<:Number}) where {S<:AbstractShape}
    _init_point_evaluation!(proxy, S)
    _set_point_evaluation!(proxy, S, shapecoordinates, ξ)
end

function _init_point_evaluation!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S}) where {Dim, RefDim, S<:AbstractShape{RefDim}}
    if proxy._shapetype[] == nnodes(S) && proxy._mode[] == 21
        return proxy
    end

    sf = proxy.shapefunc
    proxy._shapetype[] = nnodes(S)
    proxy._mode[] = 999
    proxy._quadorder[] = 0

    resize_proxy!(proxy, 1, 1, ndofs(S, sf), ndofs(S))

    return proxy
end

# TODO: Check allocations and @btime
# Normal not evaluated because we don't know where the point is located (node, edge, interior)
function _set_point_evaluation!(proxy::AbstractScalarProxy{Dim, Dim}, ::Type{S}, shapecoordinates,
        ξ::AbstractVector{<:Number}) where {Dim, S<:AbstractShape{Dim}}

    # Only for a single point ξ to avoid resizing
    @assert nnodes(S) == length(shapecoordinates)
    @assert Dim == length(ξ)

    sf = proxy.shapefunc

    proxy._evalpoints[1][1] = ξ

    for i in 1:ndofs(S, sf)
        proxy._N[1][1][i] = shapefunc_N(S, sf, i, ξ)
        proxy._dNdξ[1][1][i] = shapefunc_∇N(S, sf, i, ξ)
    end

    for i in 1:ndofs(S)
        proxy._M[1][1][i] = shapefunc_N(S, sftype(S)(), i, ξ)
        proxy._dMdξ[1][1][i] = shapefunc_∇N(S, sftype(S)(), i, ξ)
    end

    _dNdξ, _M, _dMdξ = proxy._dNdξ[1][1], proxy._M[1][1], proxy._dMdξ[1][1]

    proxy.x[1] = _to_globalcoordinates(shapecoordinates, _M)

    J = _parameterization_jacobian(shapecoordinates, _dMdξ)
    invJ = inv(J)'

    for i in eachindex(proxy.dNdx[1])
        proxy.dNdx[1][i] = invJ * _dNdξ[i]
    end

    proxy.N .= proxy._N[1]

    proxy._mode[] = 21

    return proxy
end

# Normal is evaluated because we are on a facet in higher dimensional space
function _set_point_evaluation!(proxy::AbstractScalarProxy{Dim, RefDim}, ::Type{S}, shapecoordinates,
        ξ::AbstractVector{<:Number}) where {Dim, RefDim, S<:AbstractShape{RefDim}}

    # Only for a single point ξ to avoid resizing
    @assert nnodes(S) == length(shapecoordinates)
    @assert RefDim == length(ξ)

    sf = proxy.shapefunc

    proxy._evalpoints[1][1] = ξ

    for i in 1:ndofs(S, sf)
        proxy._N[1][1][i] = shapefunc_N(S, sf, i, ξ)
    end

    for i in 1:ndofs(S)
        proxy._M[1][1][i] = shapefunc_N(S, sftype(S)(), i, ξ)
        proxy._dMdξ[1][1][i] = shapefunc_∇N(S, sftype(S)(), i, ξ)
    end

    _M, _dMdξ = proxy._M[1][1], proxy._dMdξ[1][1]

    proxy.x[1] = _to_globalcoordinates(shapecoordinates, _M)

    J = _parameterization_jacobian(shapecoordinates, _dMdξ)
    drds_x_drdt = parameterization_normal(S, J)
    norm_drds_x_drdt = norm(drds_x_drdt)

    proxy.normals[1] = drds_x_drdt / norm_drds_x_drdt

    proxy.N .= proxy._N[1]

    proxy._mode[] = 21

    return proxy
end

############################
# ... what else?
############################

############################
# Field Interpolation
############################

function interpolate(field, q::Integer)
    proxy = field.proxies[threadid()]
    interpolate(field.u, proxy, q)
end

function interpolate_cache(field, i::Integer, q::Integer)
    proxy = field.proxies[threadid()]
    interpolate(field.u_cache[i], proxy, q)
end

# TODO refactor (interpolate(u[dofs], proxy, q)?)
function interpolate(u::Vector, proxy::AbstractProxy, q::Integer)
    # @boundscheck checkbounds(u, proxy.dofs) # TODO: check only extrema(dofs)?
    # @boundscheck checkbounds(proxy.N, q)
    dofs = proxy.dofs
    N = proxy.N[q]
    @assert length(dofs) == length(N)
    u_interp = zero(first(u))
    for i in eachindex(N) # TODO: @inbounds?
        dof = dofs[i]
        u_interp += u[dof] * N[i]
    end
    return u_interp
end

# TODO More generic dispatch
function interpolate_grad(field, q::Integer)
    proxy = field.proxies[threadid()]
    interpolate_grad(field.u, proxy, q)
end

function interpolate_grad(u::Vector{T}, proxy::AbstractProxy{Dim}, q::Integer) where {Dim, T<:Number}
    # @boundscheck checkbounds(u, proxy.dofs)
    # @boundscheck checkbounds(proxy.dNdx, q)
    dofs = proxy.dofs
    ∇N = proxy.dNdx[q]
    @assert length(dofs) == length(∇N)
    ∇u_interp = zero(SVector{Dim, T})
    for i in eachindex(∇N) # TODO: @inbounds?
        dof = dofs[i]
        ∇u_interp += u[dof] * ∇N[i]
    end
    return ∇u_interp
end

function interpolate_grad(u::Vector{<:AbstractVector{T}}, proxy::AbstractProxy{Dim}, q::Integer) where {Dim, T<:Number}
    # @boundscheck checkbounds(u, proxy.dofs)
    # @boundscheck checkbounds(proxy.dNdx, q)
    dofs = proxy.dofs
    ∇N = proxy.dNdx[q]
    @assert length(dofs) == length(∇N)
    ∇u_interp = zero(SMatrix{Dim, length(first(u)), T})
    for i in eachindex(∇N) # TODO: @inbounds?
        dof = dofs[i]
        ∇u_interp += u[dof]' .* ∇N[i]
    end
    return ∇u_interp
end

function interpolate_comp_grad(field, icomp, q::Integer)
    proxy = field.proxies[threadid()]
    interpolate_comp_grad(field.u, icomp, proxy, q)
end

function interpolate_comp_grad(u::Vector{<:AbstractVector{T}}, icomp, proxy::AbstractProxy{Dim}, q::Integer) where {Dim, T<:Number}
    # @assert 1 <= icomp <= length(first(u))
    # @boundscheck checkbounds(u, proxy.dofs)
    # @boundscheck checkbounds(proxy.dNdx, q)
    dofs = proxy.dofs
    ∇N = proxy.dNdx[q]
    @assert length(dofs) == length(∇N)
    ∇u_interp = zero(SVector{Dim, T})
    for i in eachindex(∇N) # TODO: @inbounds?
        dof = dofs[i]
        ∇u_interp += u[dof][icomp] * ∇N[i]
    end
    return ∇u_interp
end
