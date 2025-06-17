using TemperFEM
using Test

@testset verbose=true "Thermal Contact 3D" begin

##############################
# Mesh
##############################

mesh_bottom = BlockMesh(0.1, 0.05, 0.05, 10, 5,  5, Hex20, origin=[0, 0, -0.025], verbose=false)
mesh_top    = BlockMesh(0.1, 0.05, 0.05, 11, 6,  5, Tet10, origin=[0, 0,  0.025], verbose=false)
mesh_ref    = BlockMesh(0.1, 0.05,  0.1, 10, 5, 10, Hex20, origin=[0.11, 0, 0], verbose=false)

##############################
# PDE
##############################

k = 30.0
cp = 1234.0
rho = 2450.0
q_bottom = -10e4
q_right = -10e4 * 0
T0 = 20.0 + 273.15

weakform = TransientHeatConduction(k, 0.0, cp, 0.0, rho, 0.0)

pde_bottom = HeatConductionProblem(weakform, mesh_bottom)
pde_bottom[] = InitialCondition(T0)

pde_top = HeatConductionProblem(weakform, mesh_top)
pde_top[] = InitialCondition(T0)

pde_ref = HeatConductionProblem(weakform, mesh_ref)
pde_ref[] = InitialCondition(T0)

##############################
# Boundary Conditions
##############################

pde_bottom[:boundary_bottom] = HeatFluxBC(q_bottom)
pde_bottom[:boundary_right] = HeatFluxBC(q_right)
pde_top[:boundary_right] = HeatFluxBC(q_right)

pde_ref[:boundary_bottom] = HeatFluxBC(q_bottom)
pde_ref[:boundary_right] = HeatFluxBC(q_right)

bc_bottom = ThermalContactBC(mesh_bottom, pde_bottom.dofmap,
    boundarytags=[:boundary_top])

bc_top = ThermalContactBC(mesh_top, pde_top.dofmap,
    boundarytags=[:boundary_bottom])

couple!(bc_bottom, bc_top; contacttol=1e-6, h_constriction=10000)

pde_top[:boundary_bottom] = bc_top
pde_bottom[:boundary_top] = bc_bottom

##############################
# Solver
##############################

# Linear solver
nprocs = 1
linsolver_bottom = MKLPardisoFactorize(;nprocs)
linsolver_top = MKLPardisoFactorize(;nprocs)
linsolver_ref = MKLPardisoFactorize(;nprocs)

# Nonlinear solver
nlsolver = NewtonSolver(xtol=0.0, rtol=1e-3, linesearch=StaticLineSearch(), verbose=false)
nlsolver_ref = NewtonSolver(xtol=0.0, rtol=1e-3, linesearch=Backtracking(), verbose=false)

# Transient solver
timesteps = PrescribedTimeSteps(range(0, 200, step=5))
solver = TransientSolver((pde_bottom, pde_top), timesteps,
    (nlsolver,), (linsolver_bottom, linsolver_top), verbose=false,
)

solver_ref = TransientSolver(pde_ref, deepcopy(timesteps), nlsolver_ref, linsolver_ref, verbose=false)

##############################
# Output
##############################
outputvtk = true
outputdir = joinpath("output", "thermal_contact_test")

vtkwriter_bottom = VtkWriter(joinpath(outputdir, "sim_bottom"), pde_bottom, enabled=outputvtk)
vtkwriter_top = VtkWriter(joinpath(outputdir, "sim_top"), pde_top, enabled=outputvtk)
vtkwriter_ref = VtkWriter(joinpath(outputdir, "sim_ref"), pde_ref, enabled=outputvtk)

##############################
# Execution
##############################

solve!(solver; writers=(vtkwriter_bottom, vtkwriter_top))
solve!(solver_ref; writers=vtkwriter_ref)

T_min = minimum(pde_top.T)
T_max = maximum(pde_bottom.T)
T_min_ref, T_max_ref = extrema(pde_ref.T)

@test T_min ≈ T_min_ref rtol=5e-3
@test T_max ≈ T_max_ref rtol=5e-3

end