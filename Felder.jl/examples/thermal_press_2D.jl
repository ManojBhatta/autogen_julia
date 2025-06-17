#=
TemperFEM.jl - Script
Date: 2024-12-19
Author: Ken Igeta
Tested in Julia 1.11.2

2D Glass Press Simulation with heat transfer between mold
and glass.

This script simulates the press model of a glass sheet between
a lower and upper mold defined in Section 2 of
"2024-12-18_TemperFEM - Thermal Press Model.pptx".
The glass is modeled as a viscoelastic material with strucural
relaxation (Narayanaswamy model). Assuming 2D plane strain.

A penalty method is used to model the contact between the glass
and the molds, including friction.

The lower mold shape is a sigmoid function with adjustable parameters.

An emperical model (Section 1) is used to model the heat transfer
between the mold and the glass.

Slide numbers "A:XX" and "B:XX" etc. refer to slide "XX" in the
following documents:

A: "2022-12_16_TemperFEM.pptx"
B: "2023-04-27_TemperFEM - Pipe and Slit Nozzle Simulation.pptx"
C: "2023-11-07_TemperFEM - Radiative Heat Transfer Model.pptx"
D: "2024-03-19_TemperFEM - Large Deformation Viscoelasticity Model.pptx"
E: "2024-07-30_TemperFEM - Press Contact Model and Results.pptx"
F: "2024-10-28_TemperFEM - Adaptive Time Stepping.pptx"
G: "2024-12-18_TemperFEM - Thermal Press Model.pptx"

A paraview state file (.pvsm) is provided in
"output/thermal_press_2D_paraviewstate.pvsm"
=#

using TemperFEM

##############################
# Parameters
##############################

# Glass dimensions
width_glass         = 0.6                       # x [m] (slide G:20)
thickness_glass     = 3e-3                      # y [m] (slide G:22)
Δx_glass            = 0.0                       # Displacement in x-direction [m])

# Mold dimensions and shape (mirrored sigmoid)
width_mold          = 0.7                       # x [m] (slide G:61)
thickness_thick     = 0.05                      # Thickness for mold shape y [m] (slide G:61)
thickness_thin      = 0.0125                    # Thickness for mold shape y [m] (slide G:61)
steepness_mold      = 100.0                     # Steepness of sigmoid function for mold shape (negative value to flip shape) [1] (slide G:61)
width_plateau       = 0.4                       # Width of plateau in sigmoid function [m] (slide G:61)
Δx_mold             = -0.02                     # Displacement in x-direction [m]

gap_glass_upper     = 0.0                       # Remaining gap between glass top and upper mold [m] >= 0
gap_glass_lower     = 0.0                       # Initial gap between glass bottom and lower mold [m] >= 0

# Upper mold translation
Δt_down             = 5.0                       # Time to move upper mold down [s] (slide G:22)
Δt_press            = 2.0                       # Time to press glass [s]  (slide G:22)
translate_y_upper   = LinearInterpolation(      # (slide B:27)
    [0.0, Δt_down, Δt_down + Δt_press, 10.0],   # Time stamps [s]  (slide G:22)
    [0.05, 0.0, 0.0, 0.05]                      # Translation in y-direction [m] (slide G:22)
)

# Glass Mesh
hx_glass            = 0.005                     # Mesh grid size in x-direction (approximately) [m] (slide G:22)
ny_glass            = 10                        # Number of elements through glass thickness (slide G:22)

# Mold Mesh
hx_mold             = 0.01                      # Mesh grid size in x-direction (approximately) [m]
ny_mold             = 4                         # Number of elements through mold thickness

# Thermal material properties (constant)
rho_glass           = 2530.0                    # Density [kg/m^3] (slide A:15, G:23)
k_th_glass          = 1.4                       # Thermal conductivity [W/m/K] (slide A:15, G:23)
cp_glass            = 1200.0                    # Specific heat [J/kg/K] (slide A:15, G:23)

# Thermal conditions as functions of position x and time t
T_init(x, t)        = 660.0 + 273.15            # Initial glass temperature [K] (slide G:20)
# T_air(x, t)         = 20.0 + 273.15           # Ambient air temperature [K] (automatically calculated using ThermalPressBC)

# Thermal convection as functions of position x and time t
h_top(x, t)         = 20.0                      # Convection coefficient [W/m^2/K] on top and edges (slide A:15, G:9)
h_bottom(x, t)      = 20.0                      # Convection coefficient [W/m^2/K] on bottom (slide A:15, G:9)
# h_edge(x, t)        = 20.0                    # Convection coefficient [W/m^2/K] on bottom (slide A:15, G:9)
g_convection_top    = 0.2                       # Gap for natural convection on top [m] (slide G:9)
g_convection_bottom = 0.2                       # Gap for natural convection on bottom [m] (slide G:9)

# Thermal Press Properties
h_constriction      = 200.0                     # Constriction conductance [W/m^2/K] in constriction region (slide G:11)
k_gap               = 0.026                     # Thermal conductivity of air [W/m/K] (slide G:10)
s_gap               = 1.0                       # Gap scaling factor [1] (slide G:10)
M_gap               = 5e-4                      # Mean gap between contacting surfaces [m] (slide G:10)
T_upper(x, t)       = 400.0 + 273.15            # Upper mold temperature [K] (slide G:9-12)
T_lower(x, t)       = 400.0 + 273.15            # Lower mold temperature [K] (slide G:9-12)

# Radiative properties
ε_glass             = [0.1, 0.4, 0.9]           # Glass emissivity (per wavelength band) [1] (slide C:16)
ϵ_mold              = 0.9                       # Mold emissivity (gray) [1] (slide C:16)
bands               = [0, 2.75e-6, 4.5e-6, Inf] # Wavelength bands for spectral properties [m] (slide C:16)

# Viscoelastic properties
bulkmod_0           = 4.1667e10                 # Bulk modulus K0 [Pa]; E = 70 GPa, ν = 0.22 (slide A:22, G:24)
bulkmod_1           = bulkmod_0 * 0.18          # Bulk modulus K1 [Pa]  (slide A:22, G:24)
shearmod            = 2.8689e10                 # Shear modulus G0 [Pa]; E = 70 GPa, ν = 0.22  (slide A:22, G:24)
C                   = [ 0.05523,  0.08205,   0.1215,   0.2286, 0.2860, 0.22662] # Structural relaxation weights [1] (slide A:11, A:26, G:24)
λ                   = [5.965e-4, 1.077e-2,   0.1362,    1.505,  6.747,   29.63] # Structural relaxation times [s] (slide A:11, A:26, G:24)
w_G                 = [ 0.05523,  0.08205,   0.1215,   0.2286, 0.2860, 0.22662] # Shear mod. Prony series weights [1] (slide A:26, A:31, G:24)
τ_G                 = [6.658e-5, 1.197e-3, 1.514e-2,   0.1672, 0.7497, 3.292]   # Shear mod. Prony series times [s] (slide A:26, A:31, G:24)
w_K                 = [0.0222,     0.0224,   0.0286,   0.2137,  0.394, 0.3191]  # Bulk mod. Prony series weights [1] (slide A:26, A:36, G:24)
τ_K                 = [5.009e-5, 9.945e-4, 2.022e-3, 1.925e-2, 0.1199, 2.033]   # Bulk mod. Prony series times [s] (slide A:26, A:36, G:24)
T_ref               = 869.0                     # Prony series reference temperature [K] (slide A:12, G:24)
HR                  = 7.62e4                    # Shift function H/R constant [1] (slide A:12, G:24)
_x                  = 0.5                       # Shift function x constant [1] (slide A:12, G:24)
Tf_init             = 550.0 + 273.15            # [K] Initial fictive temperature if glass starts as solid (= transition temperature)
α_solid             = 9e-6                      # Thermal expansion coefficient solid [1] (slide A:9, A:39, G:24)
α_liquid            = 32e-6                     # Thermal expansion coefficient liquid [1] (slide A:9, A:39, G:24)
mode                = PlaneStrain()             # 2D assumption

# Contact parameters
μ                   = 0.4                       # Friction coefficient mold-glass [1] (slide E:11)
penalty_n           = 1e10                      # Penalty factor for contact (normal) [1] (slide E:8)
penalty_t           = 1e9                       # Penalty factor for friction (tangential) [1] (slide E:12)

# Time stepping parameters
adaptive            = true                      # Use adaptive time stepping (slide F:1-14)
t_end               = 10.0                      # End time [s]
Δt_init             = 1e-5                      # Initial time step [s] (slide F:7, F:12)
Δt_min              = 1e-8                      # Minimum time step [s] (slide F:7, F:12)
Δt_max              = 0.02                      # Maximum time step [s] (slide F:7, F:12)
t_forced            = [Δt_down, Δt_down + Δt_press] # Forced time steps (to be taken by solver) [s]
ρ_lo                = 0.0                       # Minimum relative time step decrement (0 ≤ ρ_lo < 1) [1]
ρ_hi                = 1.2                       # Maximum relative time step increment (1 < ρ_hi) [1]
Δt_prescribed       = [                         # Prescribed time step phases (Time step [s], End time [s]) if adaptive = false
    (2e-3, t_end)
]

# Solver parameters
nprocs               = Threads.nthreads()       # Number of processors for PARDISO solver
rtol_nonlinear_therm = 1e-3                     # Relative tolerance for thermal nonlinear solver [1]
rtol_nonlinear_visco = 1e-3                     # Relative tolerance for viscoelastic nonlinear solver [1]
rtol_transient_therm = 1e-3                     # Relative tolerance for thermal transient solver [1] (slide F:7)
rtol_transient_visco = 1e-3                     # Relative tolerance for viscoelastic transient solver [1]  (slide F:12)

# Output options
outputdir           = joinpath("output", "thermal_press_2D") # Output directory
outputvtk           = true                      # Write .vtu files for visualization in Paraview
framestride         = 10                        # Time steps between output frames

##########################
# Model Definition
##########################

# ------ Meshing ---------

mesh_glass = BlockMesh(width_glass, thickness_glass,
    ceil(Int, width_glass/hx_glass), ny_glass,
    Quad8, origin=[0, thickness_glass/2])

mesh_lower, mesh_upper = SigmoidPressMesh(width_mold, thickness_thick, thickness_thin,
    steepness_mold, width_plateau, thickness_glass,
    ceil(Int, width_mold/hx_mold), ny_mold,
    Quad8, origin=[0, thickness_glass/2])

translate!(mesh_glass, [Δx_glass, gap_glass_lower])
translate!(mesh_lower, [Δx_mold, 0])
translate!(mesh_upper, [Δx_mold, gap_glass_upper])

# Time dependent mesh displacement for transient solver update callback
displacement_upper = TimeDependentDisplacement(mesh_upper, dy=translate_y_upper)

# writevtk(joinpath(outputdir, "mesh_glass"), mesh_glass)
writevtk(joinpath(outputdir, "mesh_mold_lower"), mesh_lower)
writevtk(joinpath(outputdir, "mesh_mold_upper"), mesh_upper)

# ------ Thermal PDE ---------

weakform_therm = TransientHeatConduction(k_th_glass, nothing, cp_glass, nothing, rho_glass, 0.0)
pde_therm = HeatConductionProblem(weakform_therm, mesh_glass, quadorder=4, quadorder_boundary=2)
pde_therm[] = InitialCondition(T_init)
# pde_therm[:boundary_top] = ThermalConvectionBC(h=h_top, T_amb=T_air)
# pde_therm[:boundary_bottom] = ThermalConvectionBC(h=h_bottom, T_amb=T_air)
# pde_therm[[:boundary_left, :boundary_right]] = ThermalConvectionBC(h=h_edge, T_amb=T_air)

pde_therm[:boundary_top] = PressThermalBC(mesh_upper, [:boundary_bottom], T_upper,
    h_convection=h_top, g_convection=g_convection_top, T_air=nothing,
    h_constriction=h_constriction,
    k_gap=k_gap, s_gap=s_gap, M_gap=M_gap,
    emissivity_source=ϵ_mold, emissivity_target=ε_glass, bands=bands)

pde_therm[:boundary_bottom] = PressThermalBC(mesh_lower, [:boundary_top], T_lower,
    h_convection=h_bottom, g_convection=g_convection_bottom, T_air=nothing,
    h_constriction=h_constriction,
    k_gap=k_gap, s_gap=s_gap, M_gap=M_gap,
    emissivity_source=ϵ_mold, emissivity_target=ε_glass, bands=bands)

# ------ Viscoelastic PDE ------

material_visco = RelaxationGlass(bulkmod_0, bulkmod_1, shearmod,
    C, λ, w_G, τ_G, w_K, τ_K,
    T_ref, HR, _x, Tf_init         , α_solid, α_liquid)
weakform_visco = Viscoelasticity2D(mode=mode, rho_0=rho_glass, gravity=[0, -9.81], model=StVenantKirchhoffModel())
pde_visco = ViscoelasticProblem{2}(weakform_visco, mesh_glass, material_visco, quadorder=4, quadorder_boundary=2)
pde_visco[] = ZeroInitialCondition()
contact_bc1 = pde_visco[:boundary_bottom] = ContactBC(mesh_lower, [:boundary_top];
    penalty_n, penalty_t, mu=μ)
contact_bc2 = pde_visco[:boundary_top] = ContactBC(mesh_upper, [:boundary_bottom];
    penalty_n, penalty_t, mu=μ)

##############################
# Solver
##############################

# Linear solver
linsolver_therm = MKLPardisoFactorize(;nprocs)
linsolver_visco = MKLPardisoFactorize(;nprocs)

# Nonlinear solver
nlsolver_therm = NewtonSolver(xtol=0.0, rtol=rtol_nonlinear_therm, linesearch=Backtracking(α_min=1e-2), maxiterations=30)
nlsolver_visco = NewtonSolver(xtol=0.0, rtol=rtol_nonlinear_visco, linesearch=Backtracking(α_min=1e-2), maxiterations=30)

# Transient solver
timesteps = if adaptive
    AdaptiveTimeSteps(t_end; Δt_init, Δt_min, Δt_max, ρ_lo, ρ_hi,
        rtol=[rtol_transient_therm, rtol_transient_visco],
        atol=[1e-6, 1e-8],
        t_forced)
else
    PrescribedTimeSteps(Δt_prescribed)
end

solver = TransientSolver((pde_therm, pde_visco), timesteps,
    (nlsolver_therm, nlsolver_visco),
    (linsolver_therm, linsolver_visco),
    callbacks=(displacement_upper,)
)

##############################
# Output
##############################

writers = if outputvtk
    (
        VtkWriter(joinpath(outputdir, "glass"), (pde_therm, pde_visco); framestride),
        BCWriter(joinpath(outputdir, "contactbc"), pde_visco, (contact_bc1, contact_bc2); framestride),
        MeshDistanceWriter(joinpath(outputdir, "mold_distance_field"), pde_visco, mesh_lower,
            [:boundary_bottom], [:boundary_top]; framestride, finalonly=true),
        MeshDistanceWriter(joinpath(outputdir, "glass_thickness_field"), pde_visco, mesh_glass,
            [:boundary_bottom], [:boundary_top]; framestride, finalonly=true),
    )
else
    nothing
end

##########################
# Execution
##########################

solve!(solver; writers)
