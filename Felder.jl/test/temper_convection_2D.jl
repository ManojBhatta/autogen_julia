using TemperFEM
using Test

# 2D tempering simulation

@testset verbose=true "Temper Convection 2D" begin

MeshTypes = [
    Tri3,
    Tri6,
    Quad4,
    Quad8,
]

matmodels = (
    LinearModel(),
    StVenantKirchhoffModel(),
)

@testset verbose=true "Adaptive time stepping: $adaptive_timestepping" for adaptive_timestepping in (false, true)

@testset verbose=true "$matmodel" for matmodel in matmodels

@testset "$MeshType" for MeshType in MeshTypes

    ##############################
    # Parameters
    ##############################

    # Glass dimensions
    width               = 0.02              # x [m]
    thickness           = 0.004             # z [m]

    # Mesh
    nz_mesh             = 16                # Number of elements through thickness (slide 37)
    h_xy_mesh           = 0.002             # Mesh grid size in xy-plane (approximately) (slide 37)

    # Constant thermal material properties
    rho                 = 2530              # Density [kg/m^3] (slide 5, 27)
    # k(T)                = 1.4               # Thermal conductivity [W/m/K] (slide 5, 27)
    # d_k(T)              = 0.0               # Thermal conductivity derivative [W/m/K^2]
    # cp(T)               = 1200.0            # Specific heat [J/kg/K] (slide 5, 27)
    # d_cp(T)             = 0.0               # Specific heat derivative [J/kg/K]

    # Thermal material properties (as a function of temperature)
    k(T)        = 0.975 + 8.58e-4 * (T - 273.15)  # Thermal conductivity [W/m/K] (slide 5, 27)
    d_k(T)      = 8.58e-4                         # Thermal conductivity derivative [W/m/K^2]
    cp(T)       = (T < 850)  ?  893.0 + 0.4 * T  :  1433.0 + 6.5e-3 * T # Specific heat [J/kg/K] (slide 5, 27)
    d_cp(T)     = (T < 850)  ?  0.4  :  6.5e-3    # Specific heat derivative [J/kg/K]

    # Thermal conditions
    T0                  = 650.0 + 273.15    # Initial glass temperature [K] (slide 5)
    T_amb               = 20.0  + 273.15    # Ambient temperature [K] (slide 5)

    # Thermal convection as functions of position x and time t (including natural convection 20 W/m^2/K)
    h_top(x, t)         = t < 10 ? 450.0 : 450.0 # Convection coefficient [W/m^2/K] on top and edges (slide 25)
    h_bottom(x, t)      = t < 10 ? 450.0 : 450.0 # Convection coefficient [W/m^2/K] on bottom (slide 25)

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
    T_trans_annealed    = 550.0 + 273.15    # [K] Annealed glass transition temperature (= initial fictive temperature Tf if T0 below transition)
    α_solid             = 9e-6              # Thermal expansion coefficient solid [1]  (slide 6, 28)
    α_liquid            = 32e-6             # Thermal expansion coefficient liquid [1]  (slide 6, 28)
    mode                = PlaneStrain()     # 2D assumption: PlaneStrain() or PlaneThermalStrain()

    # Time stepping parameters
    t_end               = 5.0               # End time [s]
    Δt_init             = 1e-5              # Initial time step [s]
    Δt_min              = 1e-8              # Minimum time step [s]
    Δt_max              = 0.1               # Maximum time step [s]
    t_forced            = [10.0]            # Forced time steps (to be taken by solver) [s]
    ρ_lo                = 0.0               # Minimum relative time step decrement (0 ≤ ρ_lo < 1) [1]
    ρ_hi                = 1.2               # Maximum relative time step increment (1 < ρ_hi) [1]
    Δt_prescribed       = [                 # Prescribed time step phases (Time step [s], End time [s]) if adaptive = false
        (20e-3, t_end),
    ]

    # Solver parameters
    nprocs               = Threads.nthreads() # Number of processors for PARDISO solver
    rtol_nonlinear_therm = 1e-3             # Relative tolerance for thermal nonlinear solver [1]
    rtol_nonlinear_visco = 1e-3             # Relative tolerance for viscoelastic nonlinear solver [1]
    rtol_transient_therm = 1e-3             # Relative tolerance for thermal transient solver [1]
    rtol_transient_visco = 1e-3             # Relative tolerance for viscoelastic transient solver [1]

    # Output options
    outputdir           = joinpath("output", "temper_convection_2D") # Output directory
    outputvtk           = false              # Write .vtu files for visualization in Paraview
    framestride         = 2                 # Time steps between output frames

    ##########################
    # Model definition
    ##########################

    # ------------------------------------ Mesh -----------------------------------
    margin = 1 * thickness                  # Refinement margin from edge
    xvals = vcat(
        linstepseries(-width/2, -width/2 + margin, thickness/nz_mesh, h_xy_mesh),
        range(-width/2 + margin, width/2 - margin, length=ceil(Int, (width - 2 * margin) / h_xy_mesh))[2:end-1],
        linstepseries(width/2 - margin, width/2, h_xy_mesh, thickness/nz_mesh)
        )
    zvals = range(-thickness/2, thickness/2, length=nz_mesh+1)
    mesh = Mesh(xvals, zvals, MeshType, verbose=false)

    # ------------------------------ Thermal PDE -----------------------------------
    weakform_therm = TransientHeatConduction(k, d_k, cp, d_cp, rho, 0)
    pde_therm = HeatConductionProblem(weakform_therm, mesh)
    pde_therm[] = InitialCondition(T0)
    pde_therm[:boundary_top]    = ThermalConvectionBC(h=h_top, T_amb=T_amb)
    pde_therm[:boundary_bottom] = ThermalConvectionBC(h=h_bottom, T_amb=T_amb)

    # ---------------------- Viscoelastic Solid Mechanics PDE ----------------------
    material_visco = RelaxationGlass(
        bulkmod_0, bulkmod_1, shearmod, C, λ, w_G, τ_G, w_K, τ_K,
        T_ref, HR, _x, T_trans_annealed, α_solid, α_liquid)
    weakform_visco = Viscoelasticity2D(rho_0=rho, mode=mode, model=matmodel)
    pde_visco = ViscoelasticProblem{2}(weakform_visco, mesh, material_visco)
    pde_visco[] = ZeroInitialCondition()
    pde_visco[:vertex_1] = ZeroDirichletVertex{2}([true, true])
    pde_visco[:vertex_2] = ZeroDirichletVertex{2}([false, true])

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

    writers = VtkWriter(joinpath(outputdir, "sim_$MeshType"), (pde_therm, pde_visco);
        framestride, enabled=outputvtk)

    ##############################
    # Execution
    ##############################

    solve!(solver; writers)

    @test minimum(first, pde_visco.σ) ≈ -4.846e+07 rtol=3e-2
    @test maximum(first, pde_visco.σ) ≈ 1.584e+07 rtol=3e-2
end

end

end

end # testset
