using TemperFEM
using Test

# Tests single body, single heater, opaque gray body, no shade, no convection

@testset verbose=true "Radiation Test 1" begin

@testset verbose=true "Adaptive time stepping: $adaptive_timestepping" for adaptive_timestepping in (false, true)

##############################
# Parameters
##############################

# Glass geometry
width_glass         = 1.0       # x [m] (slide 47)
depth_glass         = 1.0       # y [m] (slide 47)
thickness_glass     = 0.003     # z [m] (slide 47)

# Glass mesh
nx_glass            = 10        # Number of elements through width
ny_glass            = 10        # Number of elements through depth
nz_glass            = 4         # Number of elements through thickness

# Heater geometry
width_heater        = 0.5       # x [m] (slide 47)
depth_heater        = 0.5       # y [m] (slide 47)
dist_heater         = 0.5       # Distance heater to glass [m] (slide 47)

# Heater mesh
nx_heater           = 5         # Number of elements through width
ny_heater           = 5         # Number of elements through depth

# Glass thermal material properties
rho                 = 2530      # Density [kg/m^3] (slide 48)
k_th                = 1.4       # Thermal conductivity [W/m/K] (slide 48)
cp                  = 1200.0    # Specific heat [J/kg/K] (slide 48)

# Initial conditions
T0                  = 20.0 + 273.15 # Initial glass temperature [K] (slide 47)

# Radiative conditions
T_amb               = LinearInterpolation(650.0 + 273.15) # Ambient radiation temperature as function of time t

T_heater            = LinearInterpolation(800.0 + 273.15) # Heater temperature as function of time t
ε_heater            = 0.9         # Heater emissivity [1] (slide 7, 47)

ε_glass             = 0.7         # Glass emissivity [1] (slide 7, 47)
ρ_glass             = 1 - ε_glass # Glass reflectivity (opaque) [1] (slide 7, 47)

ε_ambient           = 1.0 # Ambient emissivity [1] (slide 30)

# Time stepping parameters
t_end               = 100.0             # End time [s]
Δt_init             = 1e-2              # Initial time step [s]
Δt_min              = 1e-8              # Minimum time step [s]
Δt_max              = Inf               # Maximum time step [s]
t_forced            = []                # Forced time steps (to be taken by solver) [s]
ρ_lo                = 0.0               # Minimum relative time step decrement (0 ≤ ρ_lo < 1) [1]
ρ_hi                = 1.5               # Maximum relative time step increment (1 < ρ_hi) [1]
Δt_prescribed       = range(0, 100, step=5) # Constant time step size

# Time stepping
t                   = range(0, 100, step=5) # Constant time step size (slide 48)

# Solver parameters
nprocs              = Threads.nthreads() # Number of processors for PARDISO solver
rtol_nonlinear      = 1e-3          # Relative tolerance for nonlinear solver [1]
quadorder           = 5         # View factor quadrature order (slide 27 - 40)
meshtype_glass      = Hex20     # 3D mesh type of glass (Tet4, Tet10, Hex8, Hex20)
meshtype_heater     = Hex20     # 3D mesh type of heater (Tet4, Tet10, Hex8, Hex20)

# Output options
outputdir           = joinpath("output", "radiation_test_1") # Output directory
outputvtk           = false     # Write .vtu files for visualization of Paraview
framestride         = 1                # Time steps between output frames

##############################
# Model definition
##############################

# Mesh generation
xvals = range(-width_glass/2, 0, length=nx_glass+1)
yvals = range(-depth_glass/2, 0, length=ny_glass+1)
zvals = range(-thickness_glass, 0.0, length=nz_glass+1)
mesh = Mesh(xvals, yvals, zvals, meshtype_glass, verbose=false)

# ------------------------------ Thermal PDE -----------------------------------

weakform = TransientHeatConduction(k_th, 0.0, cp, 0.0, rho, 0.0)
pde = HeatConductionProblem(weakform, mesh)
pde[] = InitialCondition(T0)

# ---------------------------- Boundary Conditions ----------------------------

ambientrad = AmbientRadiation(T_amb, emissivity=ε_ambient)

thickness_heater = 0.2
heater1 = UniformRadiosity(width_heater, depth_heater, thickness_heater,
    nx_heater, ny_heater, 1,
    T_heater, ambientrad,
    origin=[0, 0, dist_heater + thickness_heater/2],
    meshtype=meshtype_heater, boundarytags=[:boundary_bottom],
    emissivity=ε_heater,
)

radbc_top = SurfaceRadiationBC(mesh, pde.dofmap, ambientrad,
    boundarytags=[:boundary_top],
    emissivity=ε_glass, reflectivity=ρ_glass,
)

couple!(radbc_top, heater1; quadorder=quadorder)

pde[:boundary_top] = radbc_top

##############################
# Solver
##############################

# Linear solver
linsolver = MKLPardisoIterate(;nprocs)

# Nonlinear solver
nlsolver = NewtonSolver(xtol=0.0, rtol=rtol_nonlinear, verbose=false)

# Transient solver
timesteps = PrescribedTimeSteps(t)
solver = TransientSolver(pde, timesteps, nlsolver, linsolver, verbose=false)

##############################
# Output
##############################

writers = if outputvtk
    (
        VtkWriter(joinpath(outputdir, "sim"), pde; framestride),
        BCWriter(joinpath(outputdir, "sim_bcs"), pde, (heater1,); framestride),
    )
else
    nothing
end

##############################
# Execution
##############################

solve!(solver; writers)

@test minimum(pde.T) ≈ 578.4 atol=0.1
@test maximum(pde.T) ≈ 645.8 atol=0.1
@test sum(pde.T) / length(pde.T) ≈ 605.4 atol=0.1

end

end # testset

nothing