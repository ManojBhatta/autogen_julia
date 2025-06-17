#=
This script solves the weak formulation of the second order PDE

    δ ∂u/∂t - ∇ ⋅ (ν ∇u + γ u) + β ⋅ ∇u + α u = g

in 1D, 2D and 3D where the coefficients δ, ν, γ, β, α, g are
functions of the solution field u, space coordinate x and time t,
e.g. ν = ν(u, x, t). Coefficients γ and β are vectors.

Boundary conditions can be Dirichlet, Neumann or Robin type.
Backward Euler time integration is used (BDF1).

This script uses Zygote.jl for automatic differentiation (AD) to compute
the derivatives of the coefficients w.r.t. the solution field u.

The time steps of the solution are plotted live using GLMakie and
the solution is written to VTK files for post processing.

Tested in Julia 1.10 on Windows 11.
=#

using Felder
import Felder: weakformvars, residual!, jacobian!, post_step_field_update! # Explicit import to add methods
using LinearAlgebra
using StaticArrays
import Zygote

########################
# %% Plotting Setup
########################

using GLMakie: Makie, GLMakie, Figure, Axis, Axis3
using GLMakie: Label, plot, plot!, DataAspect, @lift

GLMakie.activate!(;
    float=true, # Always float on top
    focus_on_show=false, # Focus the window when newly opened
    decorated=true # Show window decorations
)

Makie.set_theme!(Makie.theme_dark())
Makie.update_theme!(backgroundcolor="#0D1117") # GitHub Dark Theme Background

########################################
# Auto Differentiation Helper Functions
########################################
# Note that Zygote AD returns `nothing` as the derivative of zero-functions,
# see https://discourse.julialang.org/t/zygote-gradient/26715

# Allocation-free workaround because Zygote allocates a Vector when using `jacobian` to
# obtain the derivative of a SVector function w.r.t scalar number, e.g. for ∂gamma/∂u and ∂beta/∂u
function svector_du_jacobian(func, u::Number, x::StaticVector{1}, t)
    f1, (df_du1,) = Zygote.withgradient((u, x, t) -> func(u, x, t)[1], u, x, t)
    return SA[f1], SA[df_du1]
end

function svector_du_jacobian(func, u::Number, x::StaticVector{2}, t)
    f1, (df_du1,) = Zygote.withgradient((u, x, t) -> func(u, x, t)[1], u, x, t)
    f2, (df_du2,) = Zygote.withgradient((u, x, t) -> func(u, x, t)[2], u, x, t)
    return SA[f1, f2], SA[df_du1, df_du2]
end

function svector_du_jacobian(func, u::Number, x::StaticVector{3}, t)
    f1, (df_du1,) = Zygote.withgradient((u, x, t) -> func(u, x, t)[1], u, x, t)
    f2, (df_du2,) = Zygote.withgradient((u, x, t) -> func(u, x, t)[2], u, x, t)
    f3, (df_du3,) = Zygote.withgradient((u, x, t) -> func(u, x, t)[3], u, x, t)
    return SA[f1, f2, f3], SA[df_du1, df_du2, df_du3]
end

########################
# %% Physics Definition
########################

@kwdef struct MyADCoefficientFormPDE{F, T1, T2, T3, T4, T5, T6} <: AbstractPDE{1, Float64}
    field::F
    delta::T1
    nu::T2
    gamma::T3
    beta::T4
    alpha::T5
    g::T6
end

function weakformvars(pde::MyADCoefficientFormPDE, fields, proxy, q)
    u = interpolate(pde.field, q)
    u_old = interpolate_cache(pde.field, 1, q)
    ∇u = interpolate_grad(pde.field, q)
    x = proxy.x[q]
    t = pde.field.t[]
    Δt = t - pde.field.t_cache[1]

    dudt = (u - u_old) / Δt # BDF1
    ddudt_du = 1 / Δt

    delta, (ddelta_du,) = Zygote.withgradient(pde.delta, u, x, t)
    nu, (dnu_du,) = Zygote.withgradient(pde.nu, u, x, t)
    gamma, dgamma_du = svector_du_jacobian(pde.gamma, u, x, t)
    beta, dbeta_du = svector_du_jacobian(pde.beta, u, x, t)
    alpha, (dalpha_du,) = Zygote.withgradient(pde.alpha, u, x, t)
    g, (dg_du,) = Zygote.withgradient(pde.g, u, x, t)

    gamma_nl = gamma
    if !isnothing(dnu_du)
        gamma_nl += dnu_du * ∇u
    end
    if !isnothing(first(dgamma_du))
        gamma_nl += dgamma_du * u
    end

    alpha_nl = alpha
    if !isnothing(first(dbeta_du))
        alpha_nl += dbeta_du ⋅ ∇u
    end
    if !isnothing(dalpha_du)
        alpha_nl += dalpha_du * u
    end
    if !isnothing(dg_du)
        alpha_nl -= dg_du
    end

    delta_nl = delta * ddudt_du
    if !isnothing(ddelta_du)
        delta_nl += ddelta_du * dudt
    end

    return (
        u = u,
        ∇u = ∇u,
        dudt = dudt,
        delta = delta,
        nu = nu,
        gamma = gamma,
        beta = beta,
        alpha = alpha,
        g = g,
        gamma_nl = gamma_nl,
        alpha_nl = alpha_nl,
        delta_nl = delta_nl,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(r, ::MyADCoefficientFormPDE, vars, proxy, i)
    @unpack u, ∇u, dudt, N, ∇N = vars
    @unpack delta, nu, gamma, beta, alpha, g = vars

    r[1] = (nu * ∇u + gamma * u) ⋅ ∇N[i] + (beta ⋅ ∇u + alpha * u - g + delta * dudt) * N[i]
end

@inline function jacobian!(J, ::MyADCoefficientFormPDE, vars, proxy, i, j)
    @unpack u, ∇u, N, ∇N = vars
    @unpack nu, gamma_nl, beta, alpha_nl, delta_nl= vars

    J[1, 1] = (nu * ∇N[j] + gamma_nl * N[j]) ⋅ ∇N[i] + (beta ⋅ ∇N[j] + (alpha_nl + delta_nl) * N[j]) * N[i]
end

function post_step_field_update!(::MyADCoefficientFormPDE, fields)
    fields.u.u_cache[1] .= fields.u.u
    fields.u.t_cache[1] = fields.u.t[]
    Makie.notify(field_node)
    sleep(0.04) # For animation
end

########################
# %% Mesh
########################

h = 0.1
mesh = Mesh(0:h:2, Edge2)
# mesh = Mesh(0:h:2, 0:h:1, Quad4) # Tri3, Tri6, Quad4, Quad8
# mesh = Mesh(0:h:2, 0:h:1, 0:h:0.5, Hex8) # Tet4, Tet10, Hex8, Hex20, Wedge6, Wedge15
# mesh = Mesh("path/to/mesh.unv")

# 2D coordinate transformation
# map!(x -> SA[cos(π/6) -sin(π/6); sin(π/6) cos(π/6)] * x, mesh.coordinates, mesh.coordinates) # Rotation
# map!(x -> SA[x[1], x[2] + 0.5 * sin(π * x[1])], mesh.coordinates, mesh.coordinates)

# 3D coordinate transformation
# map!(x -> SA[cos(π/6) -sin(π/6) 0; sin(π/6) cos(π/6) 0; 0 0 1] * x, mesh.coordinates, mesh.coordinates) # Rotation
# map!(x -> SA[x[1], x[2], x[3] + 0.5 * sin(π * x[1])], mesh.coordinates, mesh.coordinates)

########################
# %% FEM Model
########################

# PDE Coefficients: δ ∂u/∂t - ∇ ⋅ (ν ∇u + γ u) + β ⋅ ∇u + α u = g
delta(u, x, t) = 1.0
nu(u, x, t) = 0.1
gamma(u, x, t) = SA[-0.5 * t, 0.0, 0.0]
beta(u, x, t) = SA[0.0, 0.0, 0.0]
alpha(u, x, t) = exp(u) * x[1]
g(u, x, t) = x[1] < 1 ? 1.0 : 0.0

timesteps = range(0, 1, step=0.01)

fields = (
    u = ScalarField{Float64}(mesh, Lagrange(2), cachelength=1, outputname="Solution field u"),
)

pdes = (
    MyADCoefficientFormPDE(fields.u; delta, nu, gamma, beta, alpha, g),
)

bcs = (
    DirichletBC(fields.u, :boundary_left, 0.0),
    # NeumannBC(fields.u, :boundary_right, (x, t) -> -0.1 * t),
    RobinBC(fields.u, :boundary_right, (x, t) -> max(-t, -0.2), 0.0),
)

fieldhandler = FieldHandler!(fields)
bchandler = BCHandler(bcs)
soe = preassemble(pdes, bchandler)

fields.u.t_cache .= first(timesteps)

########################
# %% Plotting
########################
# Note: Makie field plots only support first order interpolation in 2D and 3D for now

fig = Figure(size=(1200, 800))
colormap = Makie.cgrad(:viridis, 17, rev=false, categorical=true) # Discrete colormap

field_node = Makie.Observable(fields.u)
axistitle_node = @lift string("Time t = ", $field_node.t[])
Makie.on(x -> Makie.autolimits!(ax1), field_node)

if ndims(mesh) == 1
    ax1 = Axis(fig[1, 1], title=axistitle_node)
    pl1 = plot!(ax1, field_node, pointsperelement=30)
    Makie.axislegend(ax1, position=:lt)
elseif ndims(mesh) == 2
    ax1 = Axis(fig[1, 1], title=axistitle_node, aspect=DataAspect())
    pl1 = plot!(ax1, field_node, colormap=colormap, strokewidth=1.0)
    Makie.Colorbar(fig[1, 2], pl1)
else
    ax1 = Axis3(fig[1, 1], title=axistitle_node, aspect=:data, perspectiveness=0.5)
    pl1 = plot!(ax1, field_node, colormap=colormap)
    Makie.Colorbar(fig[1, 2], pl1)
end

display(fig)

########################
# %% Solving
########################

writers = (
    VTKWriter("output/coefficient_form_ad_example/result", pdes, fieldhandler, order=2),
)

linsolver = LinearSolver(soe, MKLPardisoIterate(), verbose=false) # or UMFPACKFactorization, KrylovJL_BICGSTAB, KrylovJL_CG, KrylovJL_CRAIGMR, KrylovJL_GMRES, KrylovJL_LSMR, KrylovJL_MINRES, MKLPardisoFactorize, MKLPardisoIterate; see https://docs.sciml.ai/LinearSolve/stable/
nlsolver = NewtonSolver(soe, verbose=true, linsolver=linsolver, damping=1.0)
solver = TransientSolver(nlsolver, timesteps, writers=writers)

solve!(solver, soe, pdes, fieldhandler, bchandler, GaussQuadrature(5))

display(soe)
