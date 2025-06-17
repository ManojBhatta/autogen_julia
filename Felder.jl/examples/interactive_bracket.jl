using BracketSurrogateDemo
using LinearAlgebra
using StaticArrays
using Felder
using GLMakie: GLMakie, Makie, Figure, Axis, Axis3, plot, plot!, scatter!, Colorbar
using GLMakie: Observables, Observable, Label, Menu, Toggle, Button, Textbox, GridLayout
using Printf
using Cthulhu

########################
# %% GUI
########################

GLMakie.activate!(;
    inline=false, # Show plots inline in the notebook
    float=true, # Always float on top
    focus_on_show=false, # Focus the window when newly opened
    decorated=true # Show window decorations
)

# Makie.set_theme!(Makie.theme_dark())
# Makie.update_theme!(backgroundcolor="#0D1117") # GitHub Dark Theme Background

# --------------------------------

fig = Makie.Figure(size=(1280, 600), fontsize=14)

grid_settings = fig[1:2, 1] = GridLayout(width=350, height=600)

# --------------------------------

isbusy = Ref(false)

Label(grid_settings[1, 1:2], "Simulation Parameters\n", font=:bold)

Label(grid_settings[2, 1], "Young's modulus [Pa]:", halign=:left)
tb_youngsmod = Textbox(grid_settings[2, 2], stored_string="0.0", validator=Float64, halign=:left)

Label(grid_settings[3, 1], "Poisson's ratio [1]:", halign=:left)
tb_poisson = Textbox(grid_settings[3, 2], stored_string="0.0", validator=Float64, halign=:left)

Label(grid_settings[4:6, 1], "Tip traction [Pa]:", halign=:left)
tb_tractionx = Textbox(grid_settings[4, 2], stored_string="0.0", validator=Float64, halign=:left)
tb_tractiony = Textbox(grid_settings[5, 2], stored_string="0.0", validator=Float64, halign=:left)
tb_tractionz = Textbox(grid_settings[6, 2], stored_string="0.0", validator=Float64, halign=:left)

Label(grid_settings[7, 1], "Radius [m]:", halign=:left)
tb_radius = Textbox(grid_settings[7, 2], stored_string="0.0", validator=Float64, halign=:left)

Label(grid_settings[8, 1], "Mesh warp scaling:", halign=:left)
tb_warp = Textbox(grid_settings[8, 2], stored_string="1.0", validator=Float64, halign=:left)

plot1_options = [
    ("x", 1),
    ("y", 2),
    ("z", 3),
]

Label(grid_settings[9, 1], "Plot displacement comp:", halign=:left)
menu_comp = Menu(grid_settings[9, 2], options=plot1_options, default=3, width=50, halign=:left)

# --------------------------------

ax1 = Axis3(fig[2, 2],
    aspect=:data, perspectiveness=0.5, azimuth=-0.2π,
    xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",
)
ax2 = Axis3(fig[2, 3],
    aspect=:data, perspectiveness=0.5, azimuth=-0.2π,
    xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",
)

cmap1 = :coolwarm
cmap2 = :viridis

cbar1_label = Observable("?-Displacement [m]")
cbar2_label = Observable("Von Mises Stress [Pa]")

cbar1_limits = Observable((-1.0, 1.0))
cbar2_limits = Observable((0.0, 1.0))

cmap1_discrete = Makie.cgrad(cmap1, 17, rev=true, categorical=true)
cmap2_discrete = Makie.cgrad(cmap2, 17, rev=true, categorical=true)

########################
# Interaction
########################

evalparse(tb::Textbox) = eval(Meta.parse(tb.stored_string[]))

Makie.on(tb_youngsmod.stored_string) do stored_string
    youngsmod[] = evalparse(tb_youngsmod)
    update_plot(fields)
end

Makie.on(tb_poisson.stored_string) do stored_string
    poisson[] = evalparse(tb_poisson)
    update_plot(fields)
end

Makie.on(tb_tractionx.stored_string) do stored_string
    tip_traction.val[1] = evalparse(tb_tractionx)
    Observables.notify(tip_traction)
    update_plot(fields)
end

Makie.on(tb_tractiony.stored_string) do stored_string
    tip_traction.val[2] = evalparse(tb_tractiony)
    Observables.notify(tip_traction)
    update_plot(fields)
end

Makie.on(tb_tractionz.stored_string) do stored_string
    tip_traction.val[3] = evalparse(tb_tractionz)
    Observables.notify(tip_traction)
    update_plot(fields)
end

Makie.on(tb_radius.stored_string) do stored_string
    hole_radius[] = evalparse(tb_radius)
    update_plot(fields)
end

Makie.on(tb_warp.stored_string) do stored_string
    warp_scale[] = evalparse(tb_warp)
end

Makie.on(menu_comp.selection) do selection
    u_comp[] = selection
    update_plot(fields)
end

###########################
# %% FEM Model (Felder.jl)
###########################

function hole_radius_transformation(p, r_new, r_old, r_max, p0)
    @assert 0 < r_old < r_max
    @assert 0 < r_new < r_max
    x, y, z = p
    x0, y0, = p0
    r = sqrt((x - x0)^2 + (y - y0)^2)
    if r > r_max
        return p
    end
    a = r_new + (r - r_old) / (r_max - r_old) * (r_max - r_new)
    p_new = [x0 + a / r * (x - x0), y0 + a / r * (y - y0), z]
    return p_new
end

function set_hole_radius!(r)
    mesh.coordinates .= hole_radius_transformation.(mesh.initcoordinates, r, hole_r_old, hole_r_max, Ref(hole_p0))
    return mesh
end

function mises_interpolant(field, proxy, q)
    youngsmod = field.t_cache[1]
    poisson = field.t_cache[2]
    ∇u = interpolate_grad(field, q)
    ε = calc_small_strain_voigt(∇u)
    σ = hookeslaw(ε, youngsmod, poisson)
    σ_mises = calc_von_mises_stress(σ)
    return σ_mises
end

mesh_scaling = 0.1
hole_r_old = 0.2 * mesh_scaling
hole_r_max = 0.4 * mesh_scaling
hole_p0 = [0.6, 0.4] * mesh_scaling

youngsmod = Observable(1e9)
poisson = Observable(0.3)
tip_traction = Observable([0.0, 0.0, -1000])
hole_radius = Observable(0.02)

# mesh = Mesh(joinpath("meshes", "hole_bracket.unv"))
mesh = Mesh(joinpath("..", "models", "model_2.unv"))

mesh.coordinates .*= mesh_scaling
mesh.initcoordinates .= mesh.coordinates
set_hole_radius!(hole_radius[])

fields =(
    u = VectorField{3, Float64}(mesh, Lagrange(2), outputname="Displacement u"),
    mises = ScalarField(mesh, Lagrange(2), outputname="Von Mises Stress [Pa]"),
)

pdes = (
    LinearElasticity3D(fields.u),
)

bc_wall = DirichletBC(fields.u, :boundary_wall, SA[0.0, 0.0, 0.0], [true, true, true])
bc_tip = TractionBC(fields.u, [:boundary_tip], tip_traction[])

bcs = (bc_wall, bc_tip)

fieldhandler = FieldHandler!(fields)
bchandler = BCHandler(bcs)
soe = preassemble(pdes, bchandler)

solver = LinearSolver(soe, MKLPardisoFactorize(), verbose=false, update_A=true)
# solver = LinearSolver(soe, KrylovJL_CG(), verbose=true, update_A=true)
quadrature = GaussQuadrature(5)

# Abuse field cache for storing material parameters)
resize!(fields.u.t_cache, 2)
fields.u.t_cache[1] = youngsmod[]
fields.u.t_cache[2] = poisson[]

########################
# %% Plotting
########################

shapes, x_indices, u_indices = Felder._convert2isoparametricfaces1(fields.u)
coordinates = getcoordinates(mesh, SpatialFrame())
gbmesh_node = Makie.Observable(Felder.to_gb_mesh(coordinates[x_indices], shapes))

fields_node = Observable(fields)
u_node = Observable(zeros(length(u_indices)))
stress_node = Observable(zeros(length(u_indices)))
warp_scale = Observable(1000.0)
u_comp = Observable(3)

function update_plot(fields)
    comp = u_comp[]
    for i in eachindex(u_indices)
        u_node.val[i] = fields.u.u[u_indices[i]][comp]
        stress_node.val[i] = fields.mises.u[u_indices[i]]
    end
    Observables.notify(u_node)
    Observables.notify(stress_node)
    Observables.notify(warp_scale)

    cbar1_label[] = "$(("x", "y", "z")[comp])-Displacement [m]"

    u_max = maximum(abs, u_node.val)
    cbar1_limits[] = (-u_max, u_max)
    ticks1 = range(cbar1_limits.val..., length=3)
    cbar1.ticks[] = (ticks1, map(x -> @sprintf("%g", x), ticks1))


    stress_max = maximum(stress_node.val)
    cbar2_limits[] = (0.0, stress_max)
    ticks2 = range(cbar2_limits.val..., length=3)
    cbar2.ticks[] = (ticks2, map(x -> @sprintf("%g", x), ticks2))
end

Makie.on(warp_scale) do s
    for i in eachindex(u_indices)
        gbmesh_node.val.position[i] = coordinates[x_indices[i]] + s * fields.u.u[u_indices[i]]
    end
    Observables.notify(gbmesh_node)
end

Makie.on(u_comp) do comp
    update_plot(fields)
end

Observables.onany(update_plot, fields_node)

pl1 = Makie.poly!(ax1, gbmesh_node; color=u_node, colormap=cmap1, colorrange=cbar1_limits, strokewidth=1)
pl2 = Makie.poly!(ax2, gbmesh_node; color=stress_node, colormap=cmap2, colorrange=cbar2_limits, strokewidth=1)

cbar1 = Makie.Colorbar(fig[1, 2], pl1, label=cbar1_label, vertical=false)
cbar2 = Makie.Colorbar(fig[1, 3], pl2, label=cbar2_label, vertical=false)

########################
# %% Callbacks
########################

Makie.on(hole_radius) do radius
    set_hole_radius!(radius)
    solve!(solver, soe, pdes, fieldhandler, bchandler, quadrature)
    interpolate_nodalmean!(mises_interpolant, fields.mises, fields.u)
end

Makie.on(tip_traction) do traction
    bc_tip.traction .= traction
    solve!(solver, soe, pdes, fieldhandler, bchandler, quadrature)
    interpolate_nodalmean!(mises_interpolant, fields.mises, fields.u)
end

Makie.on(youngsmod) do E
    fields.u.t_cache[1] = E
    solve!(solver, soe, pdes, fieldhandler, bchandler, quadrature)
    interpolate_nodalmean!(mises_interpolant, fields.mises, fields.u)
end

Makie.on(poisson) do nu@descend foo(args...)
    fields.u.t_cache[2] = nu
    solve!(solver, soe, pdes, fieldhandler, bchandler, quadrature)
    interpolate_nodalmean!(mises_interpolant, fields.mises, fields.u)
end

function set_textbox_string!(tb, str)
    tb.stored_string.val = str
    tb.displayed_string[] = str
    return tb
end

function set_gui_elements()
    set_textbox_string!(tb_youngsmod, string(youngsmod[]))
    set_textbox_string!(tb_poisson, string(poisson[]))
    set_textbox_string!(tb_tractionx, string(tip_traction[][1]))
    set_textbox_string!(tb_tractiony, string(tip_traction[][2]))
    set_textbox_string!(tb_tractionz, string(tip_traction[][3]))
    set_textbox_string!(tb_radius, string(hole_radius[]))
    set_textbox_string!(tb_warp, string(warp_scale[]))
    return
end

function update_gui()
    set_gui_elements()
    update_plot(fields)
    return
end

########################
# Execution
########################

solve!(solver, soe, pdes, fieldhandler, bchandler, quadrature)
interpolate_nodalmean!(mises_interpolant, fields.mises, fields.u)

update_gui()

Makie.resize_to_layout!(fig)
wait(display(fig))
