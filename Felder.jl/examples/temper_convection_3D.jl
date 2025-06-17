#=
TemperFEM.jl - Script
Date: 2024-12-19
Author: Ken Igeta
Tested in Julia 1.11.2

This script solves the multiphyics model defined in
"2023-04-27_TemperFEM - Pipe and Slit Nozzle Simulation.pptx"
except that instead of a nozzle convection field a homogeneneous
convection field (for top and bottom separately) is assumed with
no motion.

Slide numbers in the comments below all refer to the above document
and provide an explanation of the parameter.

A predefined paraview state file (.pvsm) is provided in
"output/temper_convection_3D_paraviewstate.py"
=#

using TemperFEM

##############################
# Parameters
##############################

# Glass dimensions
width               = 0.1               # x [m] (slide 8)
depth               = 0.1               # y [m] (slide 8)
thickness           = 0.0031            # z [m] (slide 8)

# Mesh
nz_mesh             = 16                # Number of elements through thickness (slide 37)
h_xy_mesh           = 0.002             # Mesh grid size in xy-plane (approximately) (slide 37)

# Constant thermal material properties
rho                 = 2530              # Density [kg/m^3] (slide 5, 27)
k_th                = 1.4               # Thermal conductivity [W/m/K] (slide 5, 27)
d_k_th              = 0.0               # Thermal conductivity derivative [W/m/K^2]
cp                  = 1200.0            # Specific heat [J/kg/K] (slide 5, 27)
d_cp                = 0.0               # Specific heat derivative [J/kg/K]

# Thermal material properties (as a function of temperature)
# rho                 = 2530              # Density [kg/m^3] (slide 5, 27)
# k_th(T)             = 0.975 + 8.58e-4 * (T - 273.15)  # Thermal conductivity [W/m/K] (slide 5, 27)
# d_k_th(T)           = 8.58e-4                         # Thermal conductivity derivative [W/m/K^2]
# cp(T)               = (T < 850)  ?  893.0 + 0.4 * T  :  1433.0 + 6.5e-3 * T # Specific heat [J/kg/K] (slide 5, 27)
# d_cp(T)             = (T < 850)  ?  0.4  :  6.5e-3    # Specific heat derivative [J/kg/K]

# Thermal conditions as functions of position x and time t
T0(x, t)            = 665.0 + 273.15    # Initial glass temperature [K] (slide 5)
T_amb(x, t)         = 20.0 + 273.15     # Ambient temperature [K] (slide 5)

# Thermal convection as functions of position x and time t (including natural convection 20 W/m^2/K)
h_top(x, t)         = t < 10 ? 450.0 : 20.0 # Convection coefficient [W/m^2/K] on top and edges (slide 25)
h_bottom(x, t)      = t < 10 ? 450.0 : 20.0 # Convection coefficient [W/m^2/K] on bottom (slide 25)

# Viscoelastic properties
bulkmod_0           = 4.1666666666666664e10  # Bulk modulus K0 [Pa]; E = 70 GPa, ν = 0.22 (slide 6, 28)
bulkmod_1           = bulkmod_0 * 0.18       # Bulk modulus K1 [Pa]  (slide 6, 28)
shearmod            = 2.8688524590163937e10  # Shear modulus G0 [Pa]; E = 70 GPa, ν = 0.22  (slide 6, 28)
C                   = [ 0.05523,  0.08205, 0.1215, 0.2286, 0.2860, 0.22662]   # Structural relaxation weights [1]  (slide 6, 28)
λ                   = [5.965e-4, 1.077e-2, 0.1362,  1.505,  6.747,   29.63]   # Structural relaxation times [s]  (slide 6, 28)
w_G                 = [0.05523, 0.08205, 0.1215, 0.2286, 0.2860, 0.22662]     # Shear mod. Prony series weights [1]  (slide 6, 28)
τ_G                 = [6.658e-5, 1.197e-3, 1.514e-2, 0.1672, 0.7497, 3.292]   # Shear mod. Prony series times [s]  (slide 6, 28)
w_K                 = [0.0222, 0.0224, 0.0286, 0.2137, 0.394, 0.3191]         # Bulk mod. Prony series weights [1]  (slide 6, 28)
τ_K                 = [5.009e-5, 9.945e-4, 2.022e-3, 1.925e-2, 0.1199, 2.033] # Bulk mod. Prony series times [s]  (slide 6, 28)
T_ref               = 869.0             # Prony series reference temperature [K]  (slide 6, 28)
HR                  = 7.62e4            # Shift function H/R constant [1]  (slide 6, 28)
_x                  = 0.5               # Shift function x constant [1]  (slide 6, 28)
Tf_init             = 550.0 + 273.15    # [K] Initial fictive temperature if glass starts as solid (= transition temperature)
α_solid             = 9e-6              # Thermal expansion coefficient solid [1]  (slide 6, 28)
α_liquid            = 32e-6             # Thermal expansion coefficient liquid [1]  (slide 6, 28)

# Time stepping parameters
adaptive            = true              # Use adaptive time stepping
t_end               = 1800.0            # End time [s]
Δt_init             = 1e-5              # Initial time step [s]
Δt_min              = 1e-8              # Minimum time step [s]
Δt_max              = Inf               # Maximum time step [s]
t_forced            = []                # Forced time steps (to be taken by solver) [s]
ρ_lo                = 0.0               # Minimum relative time step decrement (0 ≤ ρ_lo < 1) [1]
ρ_hi                = 1.5               # Maximum relative time step increment (1 < ρ_hi) [1]
Δt_prescribed = [                       # (Time step [s], End time [s])
    (20e-3, 10.0),
    (0.1, 15.0),
    (1.0, 30.0),
    (5.0, 200.0),
    (20.0, 1800.0)
]

# Solver parameters
nprocs               = Threads.nthreads() # Number of processors for PARDISO solver
rtol_nonlinear_therm = 1e-3             # Relative tolerance for thermal nonlinear solver [1]
rtol_nonlinear_visco = 1e-3             # Relative tolerance for viscoelastic nonlinear solver [1]
rtol_transient_therm = 1e-3             # Relative tolerance for thermal transient solver [1]
rtol_transient_visco = 1e-3             # Relative tolerance for viscoelastic transient solver [1]

# Output options
outputdir           = joinpath("output", "temper_convection_3D") # Output directory
outputvtk           = true              # Write .vtu files for visualization in Paraview
framestride         = 1                 # Time steps between output frames

##########################
# Model definition
##########################

# ------------------------------------ Mesh -----------------------------------
margin = 5 * thickness                  # Refinement margin from edge
xvals = vcat(
    linstepseries(-width/2, -width/2 + margin, thickness/nz_mesh, h_xy_mesh),
    range(-width/2 + margin, width/2 - margin, length=ceil(Int, (width - 2 * margin) / h_xy_mesh))[2:end-1],
    linstepseries(width/2 - margin, width/2, h_xy_mesh, thickness/nz_mesh)
    )
yvals = vcat(
    linstepseries(-depth/2, -depth/2 + margin, thickness/nz_mesh, h_xy_mesh),
    range(-depth/2 + margin, depth/2 - margin, length=ceil(Int, (depth - 2 * margin) / h_xy_mesh))[2:end-1],
    linstepseries(depth/2 - margin, depth/2, h_xy_mesh, thickness/nz_mesh)
    )
zvals = range(-thickness/2, thickness/2, length=nz_mesh+1)
mesh = Mesh(xvals, yvals, zvals, Hex8)

# ------------------------------ Thermal PDE -----------------------------------
weakform_therm = TransientHeatConduction(k_th, d_k_th, cp, d_cp, rho, 0)
pde_therm = HeatConductionProblem(weakform_therm, mesh)
pde_therm[] = InitialCondition(T0)
pde_therm[:boundary_top]    = ThermalConvectionBC(h=h_top, T_amb=T_amb)
pde_therm[:boundary_bottom] = ThermalConvectionBC(h=h_bottom, T_amb=T_amb)
pde_therm[[:boundary_left, :boundary_right, :boundary_front, :boundary_back]] = ThermalConvectionBC(h=h_top, T_amb=T_amb)

# ---------------------- Viscoelastic Solid Mechanics PDE ----------------------
material_visco = RelaxationGlass(bulkmod_0, bulkmod_1, shearmod,
    C, λ, w_G, τ_G, w_K, τ_K,
    T_ref, HR, _x, Tf_init         , α_solid, α_liquid)
weakform_visco = Viscoelasticity3D(rho_0=rho, model=LinearModel())
pde_visco = ViscoelasticProblem{3}(weakform_visco, mesh, material_visco)
pde_visco[] = ZeroInitialCondition()
pde_visco[:vertex_1] = ZeroDirichletVertex{3}([true, true, true])
pde_visco[:vertex_2] = ZeroDirichletVertex{3}([false, true, true])
pde_visco[:vertex_4] = ZeroDirichletVertex{3}([false, false, true])

##############################
# Solver
##############################

# Linear solver
linsolver_therm = MKLPardisoIterate(;nprocs)
linsolver_visco = MKLPardisoIterate(;nprocs)

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
)

##############################
# Output
##############################

writers = VtkWriter(joinpath(outputdir, "sim"), (pde_therm, pde_visco); framestride, enabled=outputvtk)

##############################
# Execution
##############################

solve!(solver; writers)
