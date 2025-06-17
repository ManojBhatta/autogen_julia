using TemperFEM
using Test

# Tests single body, single heater, semi-transparent spectral glass, one shade, with convection

@testset verbose=true "Radiation Test 2" begin

##############################
# Parameters
##############################

# Glass geometry
width_glass         = 1.0       # x [m] (slide 54)
depth_glass         = 1.0       # y [m] (slide 54)
thickness_glass     = 0.003     # z [m] (slide 54)

# Glass mesh
nx_glass            = 10        # Number of elements through width
ny_glass            = 10        # Number of elements through depth
nz_glass            = 4         # Number of elements through thickness

# Heater geometry
width_heater        = 0.5       # x [m] (slide 54)
depth_heater        = 0.5       # y [m] (slide 54)
dist_heater         = 0.5       # Distance heater to glass [m] (slide 54)

# Heater mesh
nx_heater           = 10        # Number of elements through width
ny_heater           = 10        # Number of elements through depth

# Shade geometry
width_shade         = 0.5       # x [m] (slide 56)
depth_shade         = 0.5       # y [m] (slide 56)
dist_shade          = 0.25      # Distance shade to glass [m] (slide 56)
x_shade             = -0.25     # x-coordinate of shade center [m] (slide 56)
y_shade             = 0.0       # y-coordinate of shade center [m] (slide 56)

# Glass thermal material properties
rho                 = 2530      # Density [kg/m^3] (slide 56)
k_th                = 1.4       # Thermal conductivity [W/m/K] (slide 56)
cp                  = 1200.0    # Specific heat [J/kg/K] (slide 56)

# Initial conditions
T0                  = 20.0 + 273.15 # Initial glass temperature [K] (slide 54)

# Convective conditions
h_conv              = 20.0      # Convection coefficient on top surface [W/m^2/K] (slide 55)
T_air               = 640.0 + 273.15 # Air temperature [K] (slide 55)

# Radiative conditions
bands               = [0, 2.75e-6, 4.5e-6, Inf] # Wavelength bands for spectral properties [m] (slide 16, 54)
T_amb               = LinearInterpolation(650.0 + 273.15) # Ambient radiation temperature as function of time t

T_heater            = LinearInterpolation(800.0 + 273.15) # Heater temperature as function of time t
ε_heater            = 0.9         # Heater emissivity [1] (slide 7, 54)

ε_glass             = [0.1, 0.75, 0.9] # Glass emissivity [1] (slide 16, 54)
ρ_glass             = [0.1, 0.1, 0.1]  # Glass reflectivity [1] (slide 16, 54)

τ_shade             = 0.0       # Shade transmissivity [1] (slide 16, 34)

ε_ambient           = 1.0 # Ambient emissivity [1] (slide 30)

# Time stepping
t                   = range(0, 100, step=5) # Constant time step size (slide 56)

# Solver parameters
nprocs              = Threads.nthreads() # Number of processors for PARDISO solver
rtol_nonlinear      = 1e-3          # Relative tolerance for nonlinear solver [1]
quadorder           = 5         # View factor quadrature order (slide 27 - 40)
meshtype_glass      = Hex20     # 3D mesh type of glass (Tet4, Tet10, Hex8, Hex20)
meshtype_heater     = Hex20     # 3D mesh type of heater (Tet4, Tet10, Hex8, Hex20)

# Output options
outputdir           = joinpath("output", "radiation_test_2") # Output directory
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

ambientrad = AmbientRadiation(T_amb, emissivity=ε_ambient, bands=bands)

thickness_heater = 0.2
heater1 = UniformRadiosity(width_heater, depth_heater, thickness_heater,
    nx_heater, ny_heater, 1,
    T_heater, ambientrad,
    origin=[0, 0, dist_heater + thickness_heater/2],
    meshtype=meshtype_heater, boundarytags=[:boundary_bottom],
    emissivity=ε_heater, bands=bands,
)

thickness_shade = 0.1
shade1 = Shade(width_shade, depth_shade, thickness_shade,
    origin=[x_shade, y_shade, dist_shade + thickness_shade/2],
    boundarytags=[:boundary_bottom],
    transmissivity=τ_shade, bands=bands,
)

boundarytags = [:boundary_top, :boundary_bottom]

radbc_top = SurfaceRadiationBC(mesh, pde.dofmap, ambientrad,
    boundarytags=boundarytags,
    emissivity=ε_glass, reflectivity=ρ_glass, bands=bands,
)

couple!(radbc_top, heater1; shades=shade1, quadorder=quadorder)

pde[boundarytags] = radbc_top
pde[boundarytags] = ThermalConvectionBC(h=h_conv, T_amb=T_air)

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
        BCWriter(joinpath(outputdir, "sim_bcs"), pde, (heater1, shade1); framestride),
    )
else
    nothing
end

##############################
# Execution
##############################

solve!(solver; writers)

@test minimum(pde.T) ≈ 832.7 atol=0.1
@test maximum(pde.T) ≈ 848.2 atol=0.1
@test sum(pde.T) / length(pde.T) ≈ 837.7 atol=0.1

end # testset

nothing
