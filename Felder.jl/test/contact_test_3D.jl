using TemperFEM
using Test
using Statistics

# Test 3D large deformation viscoelastic model test with rigid body contact.

@testset verbose=true "Contact Test 3D" begin

MeshTypes = [
    Tet10,
    Hex20,
    Wedge15,
]

@testset "$MeshType" for MeshType in MeshTypes

    println("Running Contact Test 3D: $MeshType ...")

    ##############################
    # Parameters
    ##############################

    # Glass dimensions
    width_glass        = 0.1       # x [m]
    depth_glass        = 0.002      # y [m]
    thickness_glass    = 0.001      # z [m]

    # Mesh
    nx_glass           = 50         # Number of elements through width
    ny_glass           = 1          # Number of elements through depth
    nz_glass           = 12         # Number of elements through thickness

    # Mold geometry
    width_mold         = 0.11       # x [m]
    depth_mold         = 0.005      # y [m]
    thickness_mold     = 0.005       # z [m]
    dist_mold          = 0.000      # z [m]

    # Mold mesh
    nx_mold            = 25         # Number of elements through width
    ny_mold            = 2          # Number of elements through depth
    nz_mold            = 2          # Number of elements through thickness

    # Constant thermal material properties
    rho_glass           = 2530.0    # Density [kg/m^3]
    k_th_glass          = 1.4       # Thermal conductivity [W/m/K]
    cp_glass            = 1200.0    # Specific heat [J/kg/K]

    # Thermal conditions
    T0                  = 700.0 + 273.15    # Initial glass temperature [K]
    T_amb(x, t)         = t < 0.5 ? 700.0 + 273.15 : 400.0 + 273.15 # Ambient temperature [K]

    # Thermal convection as functions of position x and time t
    h_top(x, t)         = 500.0 # Convection coefficient [W/m^2/K] on top
    h_bottom(x, t)      = 500.0 # Convection coefficient [W/m^2/K] on bottom

    # Time stepping
    Δt_prescribed = [   # Prescribed time step phases (Time step [s], End time [s]) if adaptive = false
        (2e-3, 0.6),
        (5e-3, 1.0)
    ]

    # Solver parameters
    nprocs               = Threads.nthreads() # Number of processors for PARDISO solver
    rtol_nonlinear_therm = 1e-3             # Relative tolerance for thermal nonlinear solver [1]
    rtol_nonlinear_visco = 1e-3             # Relative tolerance for viscoelastic nonlinear solver [1]
    rtol_transient_therm = 1e-3             # Relative tolerance for thermal transient solver [1]
    rtol_transient_visco = 1e-3             # Relative tolerance for viscoelastic transient solver [1]

    # Output options
    outputdir           = joinpath("output", "contact_test_3D") # Output directory
    outputvtk           = true              # Write .vtu files for visualization in Paraview
    framestride         = 2                 # Time steps between output frames

    ##########################
    # Model definition
    ##########################

    # ------------------------------------ Mesh -----------------------------------
    # Glass Mesh generation
    xvals_g = range(-width_glass/2, width_glass/2, length=nx_glass+1)
    yvals_g = range(-depth_glass/2, depth_glass/2, length=ny_glass+1)
    zvals_g = range(-thickness_glass/2, thickness_glass/2, length=nz_glass+1)
    mesh_glass = Mesh(xvals_g, yvals_g, zvals_g, MeshType, verbose=false)

    # Mold Mesh generation
    x0_m, y0_m, z0_m = [0, 0, -thickness_mold/2 - dist_mold - thickness_glass/2]
    xvals_m = range(x0_m - width_mold/2,     x0_m + width_mold/2,     length=nx_mold+1)
    yvals_m = range(y0_m - depth_mold/2,     y0_m + depth_mold/2,     length=ny_mold+1)
    zvals_m = range(z0_m - thickness_mold/2, z0_m + thickness_mold/2, length=nz_mold+1)
    mesh_mold = Mesh(xvals_m, yvals_m, zvals_m, Hex20, verbose=false)

    # ---------------------------- Mold Shape Transformation -----------------------
    # These operations are performed on the mold mesh coordinates to obtain the desired
    # mold shape

    sigmoid(x) = 1 / (1 + exp(-x))

    function sigmoid_scaling(x, min_scale, a=500, x_flat=0.0, x_inflection=0.015)
        absx = abs(x)
        if absx < x_flat
            return 1.0
        else
            sigmoid(-a*(absx - x_inflection)) * (1 - min_scale) + min_scale
        end
    end

    function mold_shape_transform(x, z_bottom, min_scale)
        z_new = (x[3] - z_bottom) * sigmoid_scaling(x[1], min_scale) + z_bottom
        return SA[x[1], x[2], z_new]
    end

    z_bottom = -thickness_mold - dist_mold - thickness_glass/2
    for (i, x) in enumerate(mesh_mold.coordinates)
        mesh_mold.coordinates[i] = mold_shape_transform(x, z_bottom, 0.5)
    end
    mesh_mold.initcoordinates .= mesh_mold.coordinates

    # ------------------------------ Thermal PDE -----------------------------------
    weakform_therm = TransientHeatConduction(k_th_glass, 0.0, cp_glass, 0.0, rho_glass, 0.0)
    pde_therm = HeatConductionProblem(weakform_therm, mesh_glass)
    pde_therm[] = InitialCondition(T0)
    pde_therm[:boundary_top]    = ThermalConvectionBC(h=h_top, T_amb=T_amb)
    pde_therm[:boundary_bottom] = ThermalConvectionBC(h=h_bottom, T_amb=T_amb)

    # ---------------------- Viscoelastic Solid Mechanics PDE ----------------------
    material_visco = RelaxationGlass()
    weakform_visco = Viscoelasticity3D(rho_0=rho_glass, gravity=[0, 0, -9.81], model=StVenantKirchhoffModel())
    pde_visco = ViscoelasticProblem{3}(weakform_visco, mesh_glass, material_visco)
    pde_visco[] = ZeroInitialCondition()
    pde_visco[:boundary_bottom] = ContactBC(mesh_mold, [:boundary_top], penalty_n=1e10, penalty_t=1e9, mu=0.5)
    pde_visco[[:boundary_front, :boundary_back]] = ZeroDirichletBC{3}([false, true, false])

    ##########################
    # Solver
    ##########################

    # Linear solver
    linsolver_therm = MKLPardisoFactorize(;nprocs)
    linsolver_visco = MKLPardisoFactorize(;nprocs)

    # Nonlinear solver
    nlsolver_therm = NewtonSolver(xtol=0.0, rtol=rtol_nonlinear_therm, linesearch=Backtracking(α_min=1e-2), maxiterations=30, verbose=false)
    nlsolver_visco = NewtonSolver(xtol=0.0, rtol=rtol_nonlinear_visco, linesearch=Backtracking(α_min=1e-2), maxiterations=30, verbose=false)

    # Transient solver
    timesteps = PrescribedTimeSteps(Δt_prescribed)
    solver = TransientSolver((pde_therm, pde_visco), timesteps,
        (nlsolver_therm, nlsolver_visco),
        (linsolver_therm, linsolver_visco),
        verbose=false,
    )

    ##############################
    # Output
    ##############################

    writers = (
        VtkWriter(joinpath(outputdir, "sim_glass_$MeshType"), (pde_therm, pde_visco), framestride=framestride, enabled=outputvtk),
    )

    writevtk(joinpath(outputdir, "mold_mesh"), mesh_mold)

    ##########################
    # Execution
    ##########################

    solve!(solver; writers)

    @test minimum(last, pde_visco.mesh.coordinates) ≈ -3e-3 rtol=3e-2
    @test maximum(last, pde_visco.mesh.coordinates) ≈ 0.596e-3 rtol=5e-2

    @test minimum(first, pde_visco.mesh.coordinates) ≈ -0.0498 rtol=3e-2
    @test maximum(first, pde_visco.mesh.coordinates) ≈ 49.6e-3 rtol=3e-2

    @test mean(pde_visco.mesh.coordinates) ≈ [0, 0, -1.16e-3] rtol=3e-2

    # @test mean(first, pde_visco.σ) ≈ 1.12e7 rtol=5e-2
end

end

nothing
