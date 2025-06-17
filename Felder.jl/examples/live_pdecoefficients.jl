using Felder
using LinearAlgebra
using StaticArrays
import Felder: weakformvars, residual!, jacobian!, post_step_field_update!, evalargs, bcvars
using GLMakie: Makie, Figure, Axis, Axis3
using Printf

########################
# Physics Definition
########################

@kwdef struct CoefficientFormPDE{F, T1, T2, T3, T4, T5} <: AbstractPDE{1, Float64}
    field::F
    nu::T1
    gamma::T2   = SA[0.0, 0.0]
    beta::T3    = SA[0.0, 0.0]
    alpha::T4   = 0.0
    g::T5       = 0.0
end

function weakformvars(pde::CoefficientFormPDE, fields, proxy, q)
    u = interpolate(pde.field, q)
    u_old_1 = interpolate_cache(pde.field, 1, q)
    u_old_2 = interpolate_cache(pde.field, 2, q)
    x = proxy.x[q]
    t = pde.field.t[]
    Δt = t - pde.field.t_cache[1]

    dudt =  (3 * u - 4 * u_old_1 + u_old_2) / (2 * Δt) # Constant step size BDF2
    dudt_du = 3 / (2 * Δt)

    return (
        u = u,
        ∇u = interpolate_grad(pde.field, q),
        nu = evalargs(pde.nu, x, t),
        gamma = evalargs(pde.gamma, x, t),
        beta = evalargs(pde.beta, x, t),
        alpha = evalargs(pde.alpha, x, t),
        g = evalargs(pde.g, x, t),
        dudt = dudt,
        dudt_du = dudt_du,
        x = x,
        t = t,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(r, ::CoefficientFormPDE, vars, proxy, i)
    @unpack u, ∇u, nu, gamma, beta, alpha, dudt, g, N, ∇N = vars
    r[1] = (nu * ∇u + gamma * u)  ⋅ ∇N[i] + (beta ⋅ ∇u + alpha * u - g + dudt) * N[i]
end

@inline function jacobian!(J, ::CoefficientFormPDE, vars, proxy, i, j)
    @unpack u, ∇u, nu, gamma, beta, alpha, dudt, dudt_du, N, ∇N = vars
    J[1, 1] = (nu * ∇N[j] + gamma * N[j]) ⋅ ∇N[i] + (beta ⋅ ∇N[j] + (alpha + dudt_du) * N[j]) * N[i]
end

function post_step_field_update!(::CoefficientFormPDE, fields)
    fields.u.u_cache[2] .= fields.u.u_cache[1]
    fields.u.u_cache[1] .= fields.u.u
    fields.u.t_cache[1] = fields.u.t[]
end

struct MovableSource
    position::Vector{Float64}
    g::Base.RefValue{Float64}
    r::Float64
    factor::Base.RefValue{Float64}
end

gauss2D(x, σ) = exp(-(x[1]^2 + x[2]^2) / (2 * σ^2))

function evalargs(source::MovableSource, x, t)
    return gauss2D(x - source.position, source.r) * source.g[] * source.factor[]
end

########################
# FEM Model
########################

function init_fem(nx=100, ny=50)
    xmax = 0.2 # [m]
    ymax = 0.1 # [m]
    mesh = Mesh(range(0, xmax, length=(nx + 1)), range(0, ymax, length=(ny + 1)), Quad4)

    fields = (
        u = ScalarField(mesh, Lagrange(1), cachelength=2, outputname="u [°C]"),
    )

    g = MovableSource([0.0, 0.0], Ref(1.0), 0.002, Ref(0.0))

    pde = CoefficientFormPDE(fields.u, nu=Ref(3e-5), gamma=Ref(SA[0.0, 0]), beta=Ref(SA[0.0, 0]), alpha=Ref(0.0), g=g)

    bcs = (
        DirichletBC(fields.u, [:boundary_top, :boundary_bottom, :boundary_left, :boundary_right], 0.0),
    )

    fieldhandler = FieldHandler!(fields)
    bchandler = BCHandler(bcs)
    soe = preassemble(pde, bchandler)

    linsolver = LinearSolver(soe, MKLPardisoFactorize(), update_A=true)
    applydirichlet!(bchandler, fieldhandler)
    return mesh, linsolver, soe, pde, fieldhandler, bchandler, g
end

function solvestep!(linsolver, soe, pde, fieldhandler, bchandler, quadrature=GaussQuadrature(5))
    applydirichlet!(bchandler, fieldhandler)
    assemble!(soe, pde, fieldhandler, bchandler, UndefinedMaterial(), quadrature)
    solve!(linsolver, soe)
    update_solution_fields!(pde, soe)
    post_step_update!(pde, fieldhandler, bchandler)
    return soe
end

########################
# Solving and Plotting
########################

isabort = Ref(false)
isreset = Ref(false)

Makie.set_theme!(Makie.theme_dark())
fig = Makie.Figure(size=(1000, 625), fontsize=16)

Makie.Label(fig[1, 1:20], "FEM Solution of Coefficient Form PDE (BDF2 Integrator)",
    font=:bold, tellwidth=false)
Makie.Label(fig[2, 1:20], Makie.@L_str("\\frac{\\partial u}{\\partial t} - \\mathbf{\\nabla} \\cdot (\\nu \\mathbf{\\nabla} u + \\mathbf{\\gamma} u) + \\mathbf{\\beta} \\cdot \\mathbf{\\nabla} u + \\alpha u - g(x,y,t) = 0 \\quad \\text{with} \\quad u|_\\text{All boundaries} = 0"),
    font=:regular, tellwidth=false, fontsize=20)

# --------------------------------

Makie.Label(fig[3, 1], "Mesh grid size", tellwidth=false)
menu1 = Makie.Menu(fig[4, 1],
    options=[
        ("20 x 10 (h = 10 mm)", (20, 10)),
        ("40 x 20 (h = 5 mm)", (40, 20)),
        ("50 x 25 (h = 4 mm)", (50, 25)),
        ("80 x 40 (h = 2.5 mm)", (80, 40)),
        ("100 x 50 (h = 2 mm)", (100, 50)),
        ("200 x 100 (h = 1 mm)", (200, 100)),
        ("400 x 200 (h = 0.5 mm)", (400, 200)),
        ],
    default="100 x 50 (h = 2 mm)", fontsize=12,
    tellwidth=false, width=150,
    cell_color_active=Makie.RGBf(0.85, 0.85, 0.85),
)
Makie.on(menu1.selection) do s
    # println("Grid selected")
    isreset[] = true
end

Makie.Label(fig[3, 4], "Time step [s]", tellwidth=false)
tb1 = Makie.Textbox(fig[4, 4],
    stored_string="0.02", validator=Float64, tellwidth=false,# width=150,
)
Makie.on(tb1.stored_string) do s
    # println("Time step changed")
    isreset[] = true
end

# --------------------------------

Makie.Label(fig[3, 7], Makie.@L_str("\\nu"), tellwidth=false)
tb_nu = Makie.Textbox(fig[4, 7],
    stored_string="3e-5", validator=Float64, tellwidth=false,# width=150,
)

Makie.Label(fig[3, 9], Makie.@L_str("\\gamma_x"), tellwidth=false)
tb_gamma = Makie.Textbox(fig[4, 9],
    stored_string="0.0", validator=Float64, tellwidth=false,# width=150,
)

Makie.Label(fig[3, 11], Makie.@L_str("\\beta_x"), tellwidth=false)
tb_beta = Makie.Textbox(fig[4, 11],
    stored_string="0.0", validator=Float64, tellwidth=false,# width=150,
)

Makie.Label(fig[3, 13], Makie.@L_str("\\alpha"), tellwidth=false)
tb_alpha = Makie.Textbox(fig[4, 13],
    stored_string="0.0", validator=Float64, tellwidth=false,# width=150,
)

Makie.Label(fig[3, 15], Makie.@L_str("g"), tellwidth=false)
tb_g = Makie.Textbox(fig[4, 15],
    stored_string="1.0", validator=Float64, tellwidth=false,# width=150,
)

# --------------------------------

Makie.Label(fig[3, 18:20], "Color map", tellwidth=false)
cmap = Makie.Observable(Makie.cgrad(:inferno, 17, categorical=true))
clim = Makie.Observable((0.0, 0.1))

tb_clim1 = Makie.Textbox(fig[4, 18],
    stored_string=string(clim.val[1]), validator=Float64, tellwidth=true, #width=150,
)
Makie.on(tb_clim1.stored_string) do s
    clim[] = (parse(Float64, s), clim.val[2])
end

menu2 = Makie.Menu(fig[4, 19],
    options=[
        ("viridis", :viridis),
        ("inferno", :inferno),
        ("hot", :hot),
        ("coolwarm", :coolwarm),
        ("bluesreds", :bluesreds),
        ("redsblues", :redsblues),
        ("RdBu", :RdBu),
        ("Blues", :Blues ),
        ("Reds", :Reds ),
        ("jet", :jet),
        ("turbo", :turbo),
        ],
    default="inferno", fontsize=12,
    tellwidth=true, width=120,
    cell_color_active=Makie.RGBf(0.85, 0.85, 0.85),
)
Makie.on(menu2.selection) do s
    cmap[] = Makie.cgrad(s, 17, categorical=true)
end

tb_clim2 = Makie.Textbox(fig[4, 20],
    stored_string=string(clim.val[2]), validator=Float64, tellwidth=true, #width=150,
)
Makie.on(tb_clim2.stored_string) do s
    clim[] = (clim.val[1], parse(Float64, s))
end

# --------------------------------

title = Makie.Observable("t = ")
ax1 = Axis(fig[5, 1:19], xautolimitmargin=(0, 0), yautolimitmargin=(0, 0),
    title=title,
    xlabel="x [m]", ylabel="y [m]",
    aspect=Makie.DataAspect())

Makie.deregister_interaction!(ax1, :rectanglezoom)
Makie.deregister_interaction!(ax1, :dragpan)
Makie.deregister_interaction!(ax1, :limitreset)
Makie.deregister_interaction!(ax1, :scrollzoom)

Makie.Colorbar(fig[5, 20], colormap=cmap, limits=clim, label="u [1]")

display(fig)

########################
# Animation loop
########################

for i in 1:100000
    nx, ny = menu1.selection.val
    mesh, linsolver, soe, pde, fieldhandler, bchandler, g = init_fem(nx, ny)

    field_node = Makie.Observable(fieldhandler.fields.u)
    empty!(ax1)
    Makie.plot!(ax1, field_node,
        colormap=cmap,
        interpolate=true,
        colorrange=clim,
        strokewidth=0,
    )

    Makie.on(s -> (pde.nu[] = parse(Float64, tb_nu.stored_string[])), tb_nu.stored_string)
    Makie.on(s -> (pde.gamma[] = SA[parse(Float64, tb_gamma.stored_string[]), 0.0]), tb_gamma.stored_string)
    Makie.on(s -> (pde.beta[] = SA[parse(Float64, tb_beta.stored_string[]), 0.0]), tb_beta.stored_string)
    Makie.on(s -> (pde.alpha[] = parse(Float64, tb_alpha.stored_string[])), tb_alpha.stored_string)
    Makie.on(s -> (g.g[] = parse(Float64, tb_g.stored_string[])), tb_g.stored_string)

    Makie.notify(tb_nu.stored_string)
    Makie.notify(tb_gamma.stored_string)
    Makie.notify(tb_beta.stored_string)
    Makie.notify(tb_alpha.stored_string)
    Makie.notify(tb_g.stored_string)

    Δt = parse(Float64, tb1.stored_string[])
    isreset[] = false
    t = 0.0
    while !isabort[] && !isreset[]
        t += Δt
        title[] = "t = $((@sprintf "%g" t)) s (right-click to abort)"

        settime!(fieldhandler, t)
        solvestep!(linsolver, soe, pde, fieldhandler, bchandler)
        Makie.notify(field_node)

        if Makie.ispressed(ax1, Makie.Mouse.left)
            g.position .= Makie.mouseposition(ax1)
            g.factor[] = 1
        elseif Makie.ispressed(ax1, Makie.Mouse.right)
            isabort[] = true
        else
            g.factor[] = 0
        end
    end
    isabort[] && break
end
