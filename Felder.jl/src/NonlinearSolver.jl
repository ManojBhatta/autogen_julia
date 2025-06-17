export AbstractNonlinearSolver
export NewtonSolver
export nonlinear_update!
export nonlinear_field_update!

abstract type AbstractNonlinearSolver <: AbstractSolver end

###########################
# NewtonSolver
###########################

mutable struct NewtonSolver{L, W} <: AbstractNonlinearSolver
    xtol::Float64
    rtol::Float64
    maxiter::Int
    damping::Float64
    verbose::Bool
    name::String
    linsolver::L
    writers::W
end

function NewtonSolver(soe::SystemOfEquations;
        xtol=0.0,
        rtol=1e-8,
        maxiter=1000,
        damping=1.0,
        verbose=true,
        name="Newton Solver",
        linsolver=LinearSolver(soe, UMFPACKFactorization()),
        writers=nothing)

    _writers = isnothing(writers) ? nothing : tuple(writers...)

    NewtonSolver{typeof(linsolver), typeof(_writers)}(xtol, rtol, maxiter, damping, verbose, name, linsolver, _writers)
end

function solve!(solver::NewtonSolver, soe, pdes, fieldhandler, bchandler, quadrature=defaultquadrature(fieldhandler); material=UndefinedMaterial())
    @unpack xtol, rtol, maxiter, damping = solver
    @unpack linsolver, verbose, name, writers = solver

    if verbose
        printwidth = 60
        margin_left = floor(Int, (printwidth - length(name)) / 2)
        margin_right = ceil(Int, (printwidth - length(name)) / 2)
        verbose && println("-"^margin_left * " " * name * " " * "-"^margin_right)
        @printf("%-6s%16s%16s%14s%10s", "Iter.", "∥r∥_inf", "∥Δu∥₂", "Lin. niter", "Damp.")
    end

    for k in 1:maxiter
        verbose && @printf("\n%-6d", k)

        nonlinear_update!(pdes, fieldhandler, bchandler)

        applydirichlet!(bchandler, fieldhandler)

        # TODO: Only assemble Jacobian if necessary
        assemble!(soe, pdes, fieldhandler, bchandler, material, quadrature)

        rnorm = norm(soe.residual, Inf)
        verbose && @printf("%16.4e", rnorm)
        if rnorm <= rtol
            break
        end

        sol = solve!(linsolver, soe)

        update_solution_fields!(pdes, soe, damping)

        xnorm = norm(soe.Δu)
        stats = getsolverstats(sol)
        verbose && @printf("%16.4e%14d", xnorm, stats.niter)
        verbose && @printf("%10.3g", damping)
        if xnorm <= xtol
            break
        end
        if verbose && k == maxiter
            print("\nMaximum number of Newton iterations reached!")
        end
    end

    verbose && println()

    write!(writers, pdes, fieldhandler)

    return soe
end

##############################
# Update Dispatch
##############################

nonlinear_field_update!(::AbstractPDE, fields) = return

function nonlinear_update!(pdes::Union{AbstractPDE, Tuple{Vararg{AbstractPDE}}},
                           fieldhandler::AbstractFieldHandler,
                           ::BCHandler)
    nonlinear_field_update!(pdes, fieldhandler.fields)
end

function nonlinear_field_update!(pdes::Tuple{Vararg{AbstractPDE, 1}}, fields)
    pde1, = pdes
    nonlinear_field_update!(pde1, fields)
    return
end

function nonlinear_field_update!(pdes::Tuple{Vararg{AbstractPDE, 2}}, fields)
    pde1, pde2 = pdes
    nonlinear_field_update!(pde1, fields)
    nonlinear_field_update!(pde2, fields)
    return
end

function nonlinear_field_update!(pdes::Tuple{Vararg{AbstractPDE, 3}}, fields)
    pde1, pde2, pde3 = pdes
    nonlinear_field_update!(pde1, fields)
    nonlinear_field_update!(pde2, fields)
    nonlinear_field_update!(pde3, fields)
    return
end

function nonlinear_field_update!(pdes::Tuple{Vararg{AbstractPDE, 4}}, fields)
    pde1, pde2, pde3, pde4 = pdes
    nonlinear_field_update!(pde1, fields)
    nonlinear_field_update!(pde2, fields)
    nonlinear_field_update!(pde3, fields)
    nonlinear_field_update!(pde4, fields)
    return
end

function nonlinear_field_update!(pdes::Tuple{Vararg{AbstractPDE, 5}}, fields)
    pde1, pde2, pde3, pde4, pde5 = pdes
    nonlinear_field_update!(pde1, fields)
    nonlinear_field_update!(pde2, fields)
    nonlinear_field_update!(pde3, fields)
    nonlinear_field_update!(pde4, fields)
    nonlinear_field_update!(pde5, fields)
    return
end

function nonlinear_field_update!(pdes::Tuple{Vararg{AbstractPDE, 6}}, fields)
    pde1, pde2, pde3, pde4, pde5, pde6 = pdes
    nonlinear_field_update!(pde1, fields)
    nonlinear_field_update!(pde2, fields)
    nonlinear_field_update!(pde3, fields)
    nonlinear_field_update!(pde4, fields)
    nonlinear_field_update!(pde5, fields)
    nonlinear_field_update!(pde6, fields)
    return
end

function nonlinear_field_update!(pdes::Tuple{Vararg{AbstractPDE, 7}}, fields)
    pde1, pde2, pde3, pde4, pde5, pde6, pde7 = pdes
    nonlinear_field_update!(pde1, fields)
    nonlinear_field_update!(pde2, fields)
    nonlinear_field_update!(pde3, fields)
    nonlinear_field_update!(pde4, fields)
    nonlinear_field_update!(pde5, fields)
    nonlinear_field_update!(pde6, fields)
    nonlinear_field_update!(pde7, fields)
    return
end
