using BracketSurrogateDemo
using LinearAlgebra
using StaticArrays
using Felder
using GLMakie: GLMakie, Makie, Figure, Axis, Axis3, plot, plot!, scatter!, Colorbar
using GLMakie: Observables, Observable, Label, Menu, Toggle, Button, Textbox, GridLayout, SliderGrid
using Printf
using Surrogates
using ProgressMeter

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
lb_youngsmod = Label(grid_settings[2, 2], "x", halign=:left)

Label(grid_settings[3, 1], "Poisson's ratio [1]:", halign=:left)
lb_poisson = Label(grid_settings[3, 2], "x", halign=:left)

Label(grid_settings[4:6, 1], "Tip traction [Pa]:", halign=:left)
lb_tractionx = Label(grid_settings[4, 2], "0.0", halign=:left)
lb_tractiony = Label(grid_settings[5, 2], "0.0", halign=:left)
tb_tractionz = Textbox(grid_settings[6, 2], stored_string="0.0", validator=Float64, halign=:left)

Label(grid_settings[7, 1], "Radius [mm]:", halign=:left)
tb_radius = Textbox(grid_settings[7, 2], stored_string="0.0", validator=Float64, halign=:left)

Label(grid_settings[8, 1], "Mesh warp scaling:", halign=:left)
tb_warp = Textbox(grid_settings[8, 2], stored_string="1.0", validator=Float64, halign=:left)

# --------------------------------

ax1 = Axis3(fig[1, 2],
    aspect=:data, perspectiveness=0.5, azimuth=-0.2π,
    xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",
    title="FEM"
)
ax2 = Axis3(fig[1, 3],
    aspect=:data, perspectiveness=0.5, azimuth=-0.2π,
    xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",
    title="Surrogate"
)

cmap1 = :viridis
cbar1_label = Observable("Von Mises Stress [Pa]")
cbar1_limits = Observable((-1.0, 1.0))
cmap1_discrete = Makie.cgrad(cmap1, 17, rev=false, categorical=true)

title_ax1 = ax1.title
title_ax2 = ax2.title

########################
# Interaction
########################

evalparse(tb::Textbox) = eval(Meta.parse(tb.stored_string[]))

Makie.on(tb_tractionz.stored_string) do stored_string
    tip_traction.val[3] = evalparse(tb_tractionz)
    Observables.notify(tip_traction)
    update_plot(fields)
end

Makie.on(tb_radius.stored_string) do stored_string
    hole_radius[] = evalparse(tb_radius) * 1e-3
    update_plot(fields)
end

Makie.on(tb_warp.stored_string) do stored_string
    warp_scale[] = evalparse(tb_warp)
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

mesh = Mesh(joinpath("meshes", "hole_bracket.unv"))

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
quadrature = GaussQuadrature(5)

# Abuse field cache for storing material parameters)
resize!(fields.u.t_cache, 2)
fields.u.t_cache[1] = youngsmod[]
fields.u.t_cache[2] = poisson[]

########################
# %% Callbacks
########################

elapsed_t1 = Observable(0.0)
elapsed_t2 = Observable(0.0)

Makie.on(elapsed_t1) do t
    title_ax1[] = "FEM (Elapsed time: $(@sprintf("%.4f", t)) s)"
end

Makie.on(elapsed_t2) do t
    title_ax2[] = "Surrogate (Elapsed time: $(@sprintf("%.4f", t)) s)"
end

Makie.on(hole_radius) do radius
    set_hole_radius!(radius)
    elapsed_t1[] = @elapsed solve!(solver, soe, pdes, fieldhandler, bchandler, quadrature)
    interpolate_nodalmean!(mises_interpolant, fields.mises, fields.u)
    x_scaled = ((tip_traction.val[3], radius) .- lb) ./ (ub .- lb)
    elapsed_t2[] = @elapsed y_mises .= surrogate_mises(x_scaled)
end

Makie.on(tip_traction) do traction
    bc_tip.traction .= traction
    elapsed_t1[] = @elapsed solve!(solver, soe, pdes, fieldhandler, bchandler, quadrature)
    interpolate_nodalmean!(mises_interpolant, fields.mises, fields.u)
    x_scaled = ((traction[3], hole_radius[]) .- lb) ./ (ub .- lb)
    elapsed_t2[] = @elapsed y_mises .= surrogate_mises(x_scaled)
end

function set_textbox_string!(tb, str)
    tb.stored_string.val = str
    tb.displayed_string[] = str
    return tb
end

function set_gui_elements()
    lb_youngsmod.text[] = string(youngsmod[])
    lb_poisson.text[] = string(poisson[])
    set_textbox_string!(tb_tractionz, string(tip_traction[][3]))
    set_textbox_string!(tb_radius, string(1000 * hole_radius[]))
    set_textbox_string!(tb_warp, string(warp_scale[]))
    return
end

function update_gui()
    set_gui_elements()
    update_plot(fields)
    return
end

###########################
# %% Build Surrogate
###########################
# Surrogate input:  x = (tip_traction_z, hole_radius)
# Surrogate output: y = [nodal_mises_stress]
# Scaling is important for the surrogate to work properly
# "Distance" for radial basis functions only works well when the
# input variables have similar scale. Therefore, we scale the input
# tip traction and hole radius to [0, 1]
lb = [-1.5e3, 1e-3] # Lower bounds
ub = [-0.5e3, 31e-3] # Upper bounds
# # xs = sample(200, lb, ub, SobolSample())
# # xs = sample(200, lb, ub, GridSample())
xs = [((x1, x2) .- lb) ./ (ub .- lb) # Scaling
    for x1 in range(lb[1], ub[1], length=11) # Tip traction and mises stress have a linear relationship, so only 2 points should suffice? No, because RBF kernel is not a linear function?
        for x2 in range(lb[2], ub[2], length=16)]

ys_mises = Vector{Float64}[]
@showprogress desc="Creating training data..." for x in xs
    x_unscaled = x .* (ub .- lb) .+ lb
    tip_traction.val .= (0, 0, x_unscaled[1])
    hole_radius.val = x_unscaled[2]
    set_hole_radius!(hole_radius[])
    solve!(solver, soe, pdes, fieldhandler, bchandler, quadrature)
    interpolate_nodalmean!(mises_interpolant, fields.mises, fields.u)
    push!(ys_mises, deepcopy(fields.mises.u))
end

@time "Build surrogate" surrogate_mises = RadialBasis(xs, ys_mises, lb, ub, rad=linearRadial())

y_mises = zeros(length(fields.mises.u))

########################
# %% Plotting
########################

tip_traction.val .= (0.0, 0.0, -1000)
hole_radius[] = 0.02

shapes, x_indices, u_indices = Felder._convert2isoparametricfaces1(fields.u)
coordinates = getcoordinates(mesh, SpatialFrame())
gbmesh_node = Makie.Observable(Felder.to_gb_mesh(coordinates[x_indices], shapes))

fields_node = Observable(fields)
stress_node_fem = Observable(zeros(length(u_indices)))
stress_node_sur = Observable(zeros(length(u_indices)))
warp_scale = Observable(0.0)

function update_plot(fields)
    for i in eachindex(u_indices)
        stress_node_fem.val[i] = fields.mises.u[u_indices[i]]
        stress_node_sur.val[i] = y_mises[u_indices[i]]
    end
    Observables.notify(stress_node_fem)
    Observables.notify(stress_node_sur)
    Observables.notify(warp_scale)

    stress_max = maximum(stress_node_fem.val)
    cbar1_limits[] = (0.0, stress_max)
    ticks1 = range(cbar1_limits.val..., length=3)
    cbar1.ticks[] = (ticks1, map(x -> @sprintf("%g", x), ticks1))
end

Makie.on(warp_scale) do s
    for i in eachindex(u_indices)
        gbmesh_node.val.position[i] = coordinates[x_indices[i]] + s * fields.u.u[u_indices[i]]
    end
    Observables.notify(gbmesh_node)
end

Observables.onany(update_plot, fields_node)

pl1 = Makie.poly!(ax1, gbmesh_node; color=stress_node_fem, colormap=cmap1_discrete, colorrange=cbar1_limits, strokewidth=0)
pl2 = Makie.poly!(ax2, gbmesh_node; color=stress_node_sur, colormap=cmap1_discrete, colorrange=cbar1_limits, strokewidth=0)

cbar1 = Makie.Colorbar(fig[2, 2:3], pl1, label=cbar1_label, vertical=false)

########################
# %% GUI Execution
########################

solve!(solver, soe, pdes, fieldhandler, bchandler, quadrature)
interpolate_nodalmean!(mises_interpolant, fields.mises, fields.u)

update_gui()

Makie.resize_to_layout!(fig)
wait(display(fig))
