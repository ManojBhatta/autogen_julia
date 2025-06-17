export AbstractPDE
export npdes
export ncomps
export PoissonPDE
export VectorPoissonPDE
export FishersPDE
export ScalarField
export VectorField
export weakformvars
export jacobian!
export residual!
export ncomps
export cachelength

############################
# PDE Contribution
############################

struct PDEContribution{T}
    jacobian::Vector{Matrix{T}}
    residual::Vector{Vector{T}}
end

function preallocate(pdes)
    T = promote_type((dtype(pde.field) for pde in pdes)...)
    n = maximum(ncomps(pde) for pde in pdes)

    J = [zeros(T, n, n) for _ in 1:nthreads()]
    r = [zeros(T, n) for _ in 1:nthreads()]

    return PDEContribution{T}(J, r)
end

############################
# Abstract PDE
############################

abstract type AbstractPDE{N, T} end

npdes(::AbstractPDE) = 1
npdes(pdes::Tuple{Vararg{AbstractPDE}}) = length(pdes)
ncomps(::AbstractPDE{N}) where {N} = N

Base.length(::AbstractPDE) = 1
Base.iterate(pde::AbstractPDE) = (pde, nothing)
Base.iterate(::AbstractPDE, ::Any) = nothing

function (::Type{T})(field::AbstractField; kwargs...) where {T <: AbstractPDE}
    return T(; field, kwargs...)
end

cachelength(pde::AbstractPDE) = 0

function weakformvars(pde::AbstractPDE, fields, proxy, q)
    x = proxy.x[q]
    t = pde.field.t[]
    return (
        x = x,
        t = t,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(r, pde::AbstractPDE, vars, proxy, i)
    @unpack x, t, N, ∇N = vars
    error("`residual!(r, pde, vars, proxy, i)` not implemented for pde::$(typeof(pde))")
end

@inline function jacobian!(J, pde::AbstractPDE, vars, proxy, i, j)
    @unpack x, t, N, ∇N = vars
    error("`jacobian!(r, pde, vars, proxy, i, j)` not implemented for pde::$(typeof(pde))")
end

############################
# Poisson PDE
############################

"""
    PoissonPDE(f=0.0)
    PoissonPDE(f=(x, t) -> 0.0))

Poisson equation -∇⋅∇u = f(x, t).
"""
@kwdef struct PoissonPDE{F, T} <: AbstractPDE{1, Float64}
    field::F
    f::T                = 0.0
end

function weakformvars(pde::PoissonPDE, fields, proxy, q)
    x = proxy.x[q]
    t = pde.field.t[]
    return (
        ∇u = interpolate_grad(pde.field, q),
        f = evalargs(pde.f, x, t),
        x = x,
        t = t,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(r, pde::PoissonPDE, vars, proxy, i)
    @unpack ∇u, f, N, ∇N = vars
    r[1] = ∇u ⋅ ∇N[i] - f * N[i]
end

@inline function jacobian!(J, pde::PoissonPDE, vars, proxy, i, j)
    @unpack N, ∇N = vars
    J[1, 1] = ∇N[j] ⋅ ∇N[i]
end

############################
# Poisson Vector PDE
############################

"""
    VectorPoissonPDE(f1=0.0, f2=0.0)
    VectorPoissonPDE(f1=(x, t) -> 0.0), ...)

Poisson equations for a 2-component vector field u = [u1, u2].

-∇⋅∇u1 = f1(x, t)
-∇⋅∇u2 = f2(x, t)
"""
@kwdef struct VectorPoissonPDE{F, T1, T2} <: AbstractPDE{2, Float64}
    field::F
    f1::T1              = 0.0
    f2::T2              = 0.0
end

function weakformvars(pde::VectorPoissonPDE, fields, proxy, q)
    x = proxy.x[q]
    t = pde.field.t[]
    return (
        ∇u1 = interpolate_comp_grad(pde.field, 1, q),
        ∇u2 = interpolate_comp_grad(pde.field, 2, q),
        f1 = evalargs(pde.f1, x, t),
        f2 = evalargs(pde.f2, x, t),
        x = x,
        t = t,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(r, pde::VectorPoissonPDE, vars, proxy, i)
    @unpack ∇u1, ∇u2, f1, f2 = vars
    @unpack N, ∇N = vars
    r[1] = ∇u1 ⋅ ∇N[i] - f1 * N[i]
    r[2] = ∇u2 ⋅ ∇N[i] - f2 * N[i]
end

@inline function jacobian!(J, pde::VectorPoissonPDE, vars, proxy, i, j)
    @unpack N, ∇N = vars
    J[1, 1] = ∇N[j] ⋅ ∇N[i]
    J[2, 2] = ∇N[j] ⋅ ∇N[i]
end

############################
# Fisher's Equation
############################

"""
    FishersPDE(f=0.0)

Nonlinear (semilinear) Fisher's equation

    ∂u/∂t - D * ∂²u/∂x² = r * u * (1 - u).

TODO: Currently only for constant time step size.

https://en.wikipedia.org/wiki/Fisher%27s_equation
https://www.sciencedirect.com/science/article/pii/S0893965914000457
"""
@kwdef struct FishersPDE{F} <: AbstractPDE{1, Float64}
    field::F
    D::Float64
    r::Float64
end

function weakformvars(pde::FishersPDE, fields, proxy, q)
    x = proxy.x[q]
    t = pde.field.t[]
    u = interpolate(pde.field, q)
    u_old_1 = interpolate_cache(pde.field, 1, q)
    u_old_2 = interpolate_cache(pde.field, 2, q)
    Δt = t - pde.field.t_cache[1]

    dudt =  (3 * u - 4 * u_old_1 + u_old_2) / (2 * Δt) # Constant step size BDF2

    return (
        u = u,
        ∇u = interpolate_grad(pde.field, q),
        D = pde.D,
        r = pde.r,
        dudt = dudt,
        Δt = Δt,
        x = x,
        t = t,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(res, ::FishersPDE, vars, proxy, i)
    @unpack u, ∇u, D, r, dudt, Δt = vars
    @unpack N, ∇N = vars
    res[1] = dudt * N[i] + D * ∇u ⋅ ∇N[i] - r * u * (1 - u) * N[i]
end

@inline function jacobian!(Jac, ::FishersPDE, vars, proxy, i, j)
    @unpack u, ∇u, D, r, dudt, Δt = vars
    @unpack N, ∇N = vars
    Jac[1, 1] = 1.5 / Δt * N[j] * N[i] + D * ∇N[j] ⋅ ∇N[i] - r * (1 - 2 * u) * N[j] * N[i]
end

function post_step_field_update!(::FishersPDE, fields)
    fields.u.u_cache[2] .= fields.u.u_cache[1]
    fields.u.u_cache[1] .= fields.u.u
    fields.u.t_cache[1] = fields.u.t[]
end
