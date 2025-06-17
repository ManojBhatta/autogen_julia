using Felder
import Felder: weakformvars, residual!, jacobian! # Explicit import to add methods
using LinearAlgebra
using StaticArrays

########################
# %% Plotting Setup
########################

using GLMakie: Makie, GLMakie, Figure, Axis, Axis3
using GLMakie: Label, plot, plot!, DataAspect

GLMakie.activate!(;
    float=true, # Always float on top
    focus_on_show=false, # Focus the window when newly opened
    decorated=true # Show window decorations
)

Makie.set_theme!(Makie.theme_dark())
Makie.update_theme!(backgroundcolor="#0D1117") # GitHub Dark Theme Background

########################
# %% Physics Definition
########################

@kwdef struct MyCoefficientFormPDE{F, T1, T2} <: AbstractPDE{1, Float64}
    field::F
    g::T1
    dg_du::T2
end

function weakformvars(pde::MyCoefficientFormPDE, fields, proxy, q)
    u = interpolate(pde.field, q)
    ∇u = interpolate_grad(pde.field, q)
    x = proxy.x[q]
    t = pde.field.t[]

    nu = 0.1
    dnu_du = 0.0

    gamma = zero(x)
    dgamma_du = zero(x)

    beta = zero(x)
    dbeta_du = zero(x)

    alpha = exp(u) * x[1]
    dalpha_du = exp(u) * x[1]

    g = pde.g(u, x, t)
    dg_du = pde.dg_du(u, x, t)

    gamma_nl = dnu_du * ∇u + dgamma_du * u
    alpha_nl = dbeta_du ⋅ ∇u + dalpha_du * u - dg_du

    return (
        u = u,
        ∇u = ∇u,
        nu = nu,
        gamma = gamma,
        beta = beta,
        alpha = alpha,
        g = g,
        gamma_nl = gamma_nl,
        alpha_nl = alpha_nl,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(r, ::MyCoefficientFormPDE, vars, proxy, i)
    @unpack u, ∇u, N, ∇N = vars
    @unpack nu, gamma, beta, alpha, g = vars

    r[1] = (nu * ∇u + gamma * u) ⋅ ∇N[i] + (beta ⋅ ∇u + alpha * u - g) * N[i]
end

@inline function jacobian!(J, ::MyCoefficientFormPDE, vars, proxy, i, j)
    @unpack u, ∇u, N, ∇N = vars
    @unpack nu, gamma, beta, alpha = vars
    @unpack gamma_nl, alpha_nl = vars

    J[1, 1] = (nu * ∇N[j] + (gamma + gamma_nl) * N[j]) ⋅ ∇N[i] + (beta ⋅ ∇N[j] + (alpha + alpha_nl) * N[j]) * N[i]
end

########################
# %% Mesh
########################

h = 0.1
# mesh = Mesh(0:h:1, Edge2)
mesh = Mesh(0:h:2, 0:h:1, Quad4) # Tri3, Tri6, Quad4, Quad8
# mesh = Mesh(0:h:2, 0:h:1, 0:h:0.5, Hex8) # Tet4, Tet10, Hex8, Hex20, Wedge6, Wedge15
# mesh = Mesh("path/to/mesh.unv")

# plot(mesh) # 1D
# plot(mesh, axis=(aspect=DataAspect(),)) # 2D
# plot(mesh, axis=(type=Axis3, aspect=:data,)) # 3D

########################
# %% FEM Model
########################

fields = (
    u = ScalarField{Float64}(mesh, Lagrange(2), outputname="Solution field u"),
)

g(u, x, t) = 1.0 # or u
dg_du(u, x, t) = 0.0 # or 1.0

pdes = (
    MyCoefficientFormPDE(fields.u, g=g, dg_du=dg_du),
)

bcs = (
    DirichletBC(fields.u, :boundary_left, 0.0),
    NeumannBC(fields.u, :boundary_right, (x, t) -> -0.1 * t),
)

fieldhandler = FieldHandler!(fields)
bchandler = BCHandler(bcs)
soe = preassemble(pdes, bchandler)

writers = (
    VTKWriter("output/coefficient_form_example/result", pdes, fieldhandler, order=2),
)

linsolver = LinearSolver(soe, MKLPardisoIterate(), verbose=false) # or UMFPACKFactorization, KrylovJL_BICGSTAB, KrylovJL_CG, KrylovJL_CRAIGMR, KrylovJL_GMRES, KrylovJL_LSMR, KrylovJL_MINRES, MKLPardisoFactorize, MKLPardisoIterate; see https://docs.sciml.ai/LinearSolve/stable/
nlsolver = NewtonSolver(soe, verbose=true, linsolver=linsolver, damping=1.0)

t = [0, 1] # or range(0, 10, length=100)
solver = TransientSolver(nlsolver, t, writers=writers)

solve!(solver, soe, pdes, fieldhandler, bchandler, GaussQuadrature(5))

display(soe)

########################
# %% Plotting
########################
# Note: Makie field plots only support first order interpolation in 2D and 3D for now

colormap = Makie.cgrad(:viridis, 17, rev=false, categorical=true) # Discrete colormap

fig = Figure(size=(1200, 800))

if ndims(mesh) == 1
    ax1 = Axis(fig[1, 1], title="1D Result")
    pl1 = plot!(ax1, fields.u, pointsperelement=30)
    Makie.axislegend(ax1, position=:lt)
elseif ndims(mesh) == 2
    ax2 = Axis(fig[1, 1], title="2D Result", aspect=DataAspect())
    pl2 = plot!(ax2, fields.u, colormap=colormap, strokewidth=1.0)
    Makie.Colorbar(fig[1, 2], pl2)
else
    ax3 = Axis3(fig[1, 1], title="3D Result", aspect=:data, perspectiveness=0.5)
    pl3 = plot!(ax3, fields.u, colormap=colormap)
    Makie.Colorbar(fig[1, 2], pl3)
end

display(fig)
