using Felder
using LinearAlgebra
using StaticArrays
import Felder: weakformvars, residual!, jacobian!, post_step_field_update!, evalargs, bcvars
using GLMakie: Makie, Figure, Axis, Axis3
using Printf

########################
# Physics Definition
########################

@kwdef struct ConvectionDiffusionPDE{F, T1, T2, T3} <: AbstractPDE{1, Float64}
    field::F
    D::T1
    g::T2       = 0.0
    v::T3       = SA[0.0, 0.0]
end

function weakformvars(pde::ConvectionDiffusionPDE, fields, proxy, q)
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
        D = evalargs(pde.D, x, t),
        g = evalargs(pde.g, x, t),
        v = evalargs(pde.v, x, t),
        dudt = dudt,
        dudt_du = dudt_du,
        x = x,
        t = t,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(r, ::ConvectionDiffusionPDE, vars, proxy, i)
    @unpack u, ∇u, D, g, v, dudt, N, ∇N = vars
    r[1] = D * ∇u ⋅ ∇N[i] + (v ⋅ ∇u - g + dudt) * N[i]
end

@inline function jacobian!(J, ::ConvectionDiffusionPDE, vars, proxy, i, j)
    @unpack u, ∇u, D, v, dudt, dudt_du, N, ∇N  = vars
    J[1, 1] = D * ∇N[j] ⋅ ∇N[i] + (v ⋅ ∇N[j] + dudt_du * N[j]) * N[i]
end

function post_step_field_update!(::ConvectionDiffusionPDE, fields)
    fields.u.u_cache[2] .= fields.u.u_cache[1]
    fields.u.u_cache[1] .= fields.u.u
    fields.u.t_cache[1] = fields.u.t[]
end

#------------------------------------
# ↓ Usage unclear, results strange
# struct FreeDiffusionBC{F, T1} <: AbstractIntegratedBC
#     field::F
#     boundaries::Vector{Symbol}
#     D::T1
# end

# function FreeDiffusionBC(field, boundarytags, D)
#     _btags = boundarytags isa Symbol ? [boundarytags] : boundarytags
#     FreeDiffusionBC{typeof(field), typeof(D)}(field, _btags, D)
# end

# function bcvars(bc::FreeDiffusionBC, fields, x, t, n, material, q)
#     return (
#         ∇u = interpolate_grad(bc.field, q),
#         D = evalargs(bc.D, x, t),
#     )
# end

# @inline function residual!(r, ::FreeDiffusionBC, vars, x, t, n, N, ∇N, i)
#     @unpack ∇u, D = vars
#     r[1] = -D * (∇u ⋅ n) * N[i]
# end

# @inline function jacobian!(J, ::FreeDiffusionBC, vars, x, t, n, N, ∇N, i, j)
#     @unpack ∇u, D = vars
#     J[1, 1] = -D * (∇N[j] ⋅ n) * N[i]
# end
#------------------------------------

struct MovableSource
    position::Vector{Float64}
    q::Float64
    r::Float64
    factor::Base.RefValue{Float64}
end

gauss2D(x, σ) = exp(-(x[1]^2 + x[2]^2) / (2 * σ^2))

function evalargs(source::MovableSource, x, t)
    return gauss2D(x - source.position, source.r) * source.q * source.factor[]
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

    g = MovableSource([0.0, 0.0], 1.0, 0.002, Ref(0.0))

    pde = ConvectionDiffusionPDE(fields.u, D=Ref(3e-5), g=g, v=Ref(SA[1e-2, 0]))

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

Makie.Label(fig[1, 1:20], "FEM Solution of Convection-Diffusion Equation (BDF2 Integrator)",
    font=:bold, tellwidth=false)
Makie.Label(fig[2, 1:20], Makie.@L_str("\\frac{\\partial u}{\\partial t}  = \\mathbf{\\nabla} \\cdot (D \\mathbf{\\nabla} u) - \\mathbf{v} \\cdot \\mathbf{\\nabla} u + g(x,y,t) \\quad \\text{with} \\quad u|_\\text{All boundaries} = 0"),
    font=:regular, tellwidth=false, fontsize=20)

Makie.Label(fig[3, 1:2], "Mesh grid size", tellwidth=false)
menu1 = Makie.Menu(fig[4, 1:2],
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
    tellwidth=true, width=160,
    cell_color_active=Makie.RGBf(0.85, 0.85, 0.85),
)

Makie.on(menu1.selection) do s
    # println("Grid selected")
    isreset[] = true
end

Makie.Label(fig[3, 5], "Time step [s]", tellwidth=false)
tb1 = Makie.Textbox(fig[4, 5],
    stored_string="0.02", validator=Float64, tellwidth=false,# width=150,
)
Makie.on(tb1.stored_string) do s
    # println("Time step changed")
    isreset[] = true
end

sg1 = Makie.SliderGrid(fig[3:4, 7:18],
    (label=Makie.@L_str("Diffusivity \$D\$ [m²/s]"), range=1e-6:1e-6:1e-4, startvalue=3e-5, snap=false),
    (label=Makie.@L_str("Velocity \$v_x\$ [m/s]"), range=0:1e-3:1.0, startvalue=0, snap=false),
)

clim = Makie.Observable((0.0, 0.1))
Makie.Label(fig[3, 20], "Color limit", tellwidth=false)
tb2 = Makie.Textbox(fig[4, 20],
    stored_string=string(clim.val[2]), validator=Float64, tellwidth=false, #width=150,
)
Makie.on(tb2.stored_string) do s
    # println("Cmax changed")
    clim[] = (0.0, parse(Float64, s))
end

title = Makie.Observable("t = ")
ax1 = Axis(fig[5, 1:19], xautolimitmargin=(0, 0), yautolimitmargin=(0, 0),
    title=title,
    xlabel="x [m]", ylabel="y [m]",
    aspect=Makie.DataAspect())

Makie.deregister_interaction!(ax1, :rectanglezoom)
Makie.deregister_interaction!(ax1, :dragpan)
Makie.deregister_interaction!(ax1, :limitreset)
Makie.deregister_interaction!(ax1, :scrollzoom)

# colormap=:hot
# colormap=:viridis
colormap=:inferno
# colormap=:YlOrRd
# colormap=:Blues

cmap = Makie.cgrad(colormap, 17, rev=false, categorical=true)

Makie.Colorbar(fig[5, 20]; colormap=cmap, limits=clim, label="u [1]")

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

    Δt = parse(Float64, tb1.stored_string[])
    isreset[] = false
    t = 0.0
    while !isabort[] && !isreset[]
        t += Δt
        title[] = "t = $((@sprintf "%g" t)) s (right-click to abort)"

        pde.D[] = sg1.sliders[1].value.val
        pde.v[] = SA[sg1.sliders[2].value.val, 0.0]

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
