using TemperFEM
using Test
using Statistics

# Test 2D large deformation viscoelastic model test with rigid body contact.

@testset verbose=true "Contact Test 2D" begin

MeshTypes = [
    Tri6,
    Quad8,
]

@testset verbose=true "Adaptive time stepping: $adaptive_timestepping" for adaptive_timestepping in (false, true)
    @testset verbose=true "$MeshType" for MeshType in MeshTypes
        ##############################
        # Parameters
        ##############################

        # Glass dimensions
        width_glass        = 0.1       # x [m]
        thickness_glass    = 0.001     # z [m]

        # Mesh
        nx_glass           = 50         # Number of elements through width
        nz_glass           = 12         # Number of elements through thickness

        # Mold geometry
        width_mold         = 0.11       # x [m]
        thickness_mold     = 0.005       # z [m]
        dist_mold          = 0.000      # z [m]

        # Mold mesh
        nx_mold            = 25         # Number of elements through width
        nz_mold            = 2          # Number of elements through thickness

        # Constant thermal material properties
        rho_glass           = 2530.0    # Density [kg/m^3]
        k_th_glass          = 1.4       # Thermal conductivity [W/m/K]
        cp_glass            = 1200.0    # Specific heat [J/kg/K]

        # Thermal conditions
        T0                  = 700.0 + 273.15    # Initial glass temperature [K]
        T_amb(x, t)         = t < 0.5  ? 700.0 + 273.15 : 400.0 + 273.15     # Ambient temperature [K]

        # Thermal convection as functions of position x and time t
        h_top(x, t)         = 500.0 # Convection coefficient [W/m^2/K] on top
        h_bottom(x, t)      = 500.0 # Convection coefficient [W/m^2/K] on bottom

        # Time stepping parameters
        t_end               = 1.0               # End time [s]
        Δt_init             = 1e-5              # Initial time step [s]
        Δt_min              = 1e-8              # Minimum time step [s]
        Δt_max              = 0.1               # Maximum time step [s]
        t_forced            = [0.5]             # Forced time steps (to be taken by solver) [s]
        ρ_lo                = 0.0               # Minimum relative time step decrement (0 ≤ ρ_lo < 1) [1]
        ρ_hi                = 1.2               # Maximum relative time step increment (1 < ρ_hi) [1]
        Δt_prescribed       = [                 # Prescribed time step phases (Time step [s], End time [s]) if adaptive = false
            (2e-3, t_end),
        ]

        # Solver parameters
        nprocs               = Threads.nthreads() # Number of processors for PARDISO solver
        rtol_nonlinear_therm = 1e-3             # Relative tolerance for thermal nonlinear solver [1]
        rtol_nonlinear_visco = 1e-3             # Relative tolerance for viscoelastic nonlinear solver [1]
        rtol_transient_therm = 1e-3             # Relative tolerance for thermal transient solver [1]
        rtol_transient_visco = 1e-3             # Relative tolerance for viscoelastic transient solver [1]

        # Output options
        outputdir           = joinpath("output", "contact_test_2D") # Output directory
        outputvtk           = true              # Write .vtu files for visualization in Paraview
        framestride         = 2                 # Time steps between output frames

        ##########################
        # Model definition
        ##########################

        # ------------------------------------ Mesh -----------------------------------
        # Glass Mesh generation
        xvals_g = range(-width_glass/2, width_glass/2, length=nx_glass+1)
        yvals_g = range(-thickness_glass/2, thickness_glass/2, length=nz_glass+1)
        mesh_glass = Mesh(xvals_g, yvals_g, MeshType, verbose=false)

        # Mold Mesh generation
        x0_m, y0_m = [0, -thickness_mold/2 - dist_mold - thickness_glass/2]
        xvals_m = range(x0_m - width_mold/2,     x0_m + width_mold/2,     length=nx_mold+1)
        yvals_m = range(y0_m - thickness_mold/2,     y0_m + thickness_mold/2, length=nz_mold+1)
        mesh_mold = Mesh(xvals_m, yvals_m, Quad8, verbose=false)

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

        function mold_shape_transform(x, y_bottom, min_scale)
            y_new = (x[2] - y_bottom) * sigmoid_scaling(x[1], min_scale) + y_bottom
            return SA[x[1], y_new]
        end

        y_bottom = -thickness_mold - dist_mold - thickness_glass/2
        for (i, x) in enumerate(mesh_mold.coordinates)
            mesh_mold.coordinates[i] = mold_shape_transform(x, y_bottom, 0.5)
        end
        mesh_mold.initcoordinates .= mesh_mold.coordinates

        # ------------------------------ Thermal PDE -----------------------------------
        weakform_therm = TransientHeatConduction(k_th_glass, nothing, cp_glass, nothing, rho_glass, 0.0)
        pde_therm = HeatConductionProblem(weakform_therm, mesh_glass)
        pde_therm[] = InitialCondition(T0)
        pde_therm[:boundary_top]    = ThermalConvectionBC(h=h_top, T_amb=T_amb)
        pde_therm[:boundary_bottom] = ThermalConvectionBC(h=h_bottom, T_amb=T_amb)

        # ---------------------- Viscoelastic Solid Mechanics PDE ----------------------
        material_visco = RelaxationGlass()
        weakform_visco = Viscoelasticity2D(rho_0=rho_glass, gravity=[1.0, -9.81], model=StVenantKirchhoffModel())
        pde_visco = ViscoelasticProblem{2}(weakform_visco, mesh_glass, material_visco)
        pde_visco[] = ZeroInitialCondition()
        pde_visco[:boundary_bottom] = ContactBC(mesh_mold, [:boundary_top];
            penalty_n=1e10, penalty_t=1e9, mu=0.5,
        )

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
        timesteps = if adaptive_timestepping
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
            verbose=false,
        )

        ##############################
        # Output
        ##############################

        writers = (
            VtkWriter(joinpath(outputdir, "sim_glass_$MeshType"), (pde_therm, pde_visco); framestride, enabled=outputvtk),
        )

        writevtk(joinpath(outputdir, "mold_mesh"), mesh_mold)

        ##########################
        # Execution
        ##########################

        solve!(solver; writers)

        @test minimum(last, pde_visco.mesh.coordinates) ≈ -3e-3 rtol=3e-2
        @test maximum(last, pde_visco.mesh.coordinates) ≈ 0.596e-3 rtol=3e-2

        @test minimum(first, pde_visco.mesh.coordinates) ≈ -0.0498 rtol=3e-2
        @test maximum(first, pde_visco.mesh.coordinates) ≈ 49.6e-3 rtol=3e-2

        @test mean(pde_visco.mesh.coordinates) ≈ [2.036e-5, -1.16e-3] rtol=3e-2

        # @test mean(first, pde_visco.σ) ≈ -374.4e3 rtol=3e-2
    end
end

end

nothing