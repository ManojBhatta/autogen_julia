export AbstractTransientSolver
export TransientSolver
export SteadyStateSolver
export pre_step_update!
export post_step_update!
export pre_step_field_update!
export post_step_field_update!

abstract type AbstractTransientSolver <: AbstractSolver end

###########################
# TransientSolver
###########################
# TODO: Implement variable time step schemes and strategies

mutable struct TransientSolver{T, W} <: AbstractTransientSolver
    stepsolver::T
    timesteps::Vector{Float64}
    verbose::Bool
    name::String

    writers::W
    framestride::Int
    _stridecounter::Int
end

function TransientSolver(stepsolver::T, timesteps;
        verbose=true,
        name="Transient Solver",
        writers=nothing, nframes=300,
        ) where T<:AbstractSolver

    _timesteps = collect(timesteps)

    _writers = isnothing(writers) ? nothing : tuple(writers...)
    framestride = max(1, length(_timesteps) ÷ nframes)
    _stridecounter = framestride

    @assert framestride > 0 "`framestride` must be greater than zero, e.g. 1 for every
        frame, 2 for every second frame ..."

    TransientSolver{T, typeof(_writers)}(
        stepsolver,
        _timesteps,
        verbose,
        name,
        _writers,
        framestride,
        _stridecounter,
    )
end

function TransientSolver(::SystemOfEquations, stepsolver::AbstractSolver, timesteps; kwargs...)
    TransientSolver(stepsolver, timesteps; kwargs...)
end

function TransientSolver(soe::SystemOfEquations, timesteps; verbose=true, kwargs...)
    defaultsolver = NewtonSolver(soe; linsolver=LinearSolver(soe, UMFPACKFactorization()), verbose)
    TransientSolver(defaultsolver, timesteps; kwargs...)
end

function SteadyStateSolver(stepsolver::AbstractSolver; kwargs...)
    TransientSolver(stepsolver, [0.0, Inf]; name="Steady-State Solver", kwargs...)
end

function SteadyStateSolver(soe::SystemOfEquations; verbose=true, kwargs...)
    defaultsolver = NewtonSolver(soe; linsolver=LinearSolver(soe, UMFPACKFactorization()), verbose)
    SteadyStateSolver(defaultsolver; kwargs...)
end

function solve!(solver::TransientSolver, soe, pdes, fieldhandler, bchandler, quadrature=defaultquadrature(fieldhandler);
                material=UndefinedMaterial())
    @unpack stepsolver, timesteps, verbose, name, writers = solver
    nsteps = length(timesteps) - 1
    hline = "="^62

    write!(solver, writers, pdes, fieldhandler, first(timesteps))

    for i in 1:nsteps
        t = timesteps[i + 1]
        verbose && println(hline)
        verbose && @printf("Timestep: %8d / %-8d  |  t = %.16g\n", i, nsteps, t)

        settime!(fieldhandler.fields, t)
        pre_step_update!(pdes, fieldhandler, bchandler)

        # TODO: Check if step converged
        solve!(stepsolver, soe, pdes, fieldhandler, bchandler, quadrature; material)

        write!(solver, writers, pdes, fieldhandler, t)

        post_step_update!(pdes, fieldhandler, bchandler)
    end

    if solver._stridecounter > 0
        write!(writers, pdes, fieldhandler, last(timesteps))
    end

    return soe
end

function write!(solver::TransientSolver, writers, args...; kwargs...)
    solver._stridecounter += 1
    if solver._stridecounter >= solver.framestride
        solver._stridecounter = 0
        write!(writers, args...; kwargs...)
    end
    return
end

##############################
# Update Dispatch
##############################

pre_step_field_update!(::AbstractPDE, fields) = return

function pre_step_update!(pdes, fieldhandler::AbstractFieldHandler, ::BCHandler)
    pre_step_field_update!(pdes, fieldhandler.fields)
end

function pre_step_field_update!(pdes::TupleOfLength{1, AbstractPDE}, fields)
    pde1, = pdes
    pre_step_field_update!(pde1, fields)
    return
end

function pre_step_field_update!(pdes::TupleOfLength{2, AbstractPDE}, fields)
    pde1, pde2 = pdes
    pre_step_field_update!(pde1, fields)
    pre_step_field_update!(pde2, fields)
    return
end

function pre_step_field_update!(pdes::TupleOfLength{3, AbstractPDE}, fields)
    pde1, pde2, pde3 = pdes
    pre_step_field_update!(pde1, fields)
    pre_step_field_update!(pde2, fields)
    pre_step_field_update!(pde3, fields)
    return
end

function pre_step_field_update!(pdes::TupleOfLength{4, AbstractPDE}, fields)
    pde1, pde2, pde3, pde4 = pdes
    pre_step_field_update!(pde1, fields)
    pre_step_field_update!(pde2, fields)
    pre_step_field_update!(pde3, fields)
    pre_step_field_update!(pde4, fields)
    return
end

function pre_step_field_update!(pdes::TupleOfLength{5, AbstractPDE}, fields)
    pde1, pde2, pde3, pde4, pde5 = pdes
    pre_step_field_update!(pde1, fields)
    pre_step_field_update!(pde2, fields)
    pre_step_field_update!(pde3, fields)
    pre_step_field_update!(pde4, fields)
    pre_step_field_update!(pde5, fields)
    return
end

function pre_step_field_update!(pdes::TupleOfLength{6, AbstractPDE}, fields)
    pde1, pde2, pde3, pde4, pde5, pde6 = pdes
    pre_step_field_update!(pde1, fields)
    pre_step_field_update!(pde2, fields)
    pre_step_field_update!(pde3, fields)
    pre_step_field_update!(pde4, fields)
    pre_step_field_update!(pde5, fields)
    pre_step_field_update!(pde6, fields)
    return
end

function pre_step_field_update!(pdes::TupleOfLength{7, AbstractPDE}, fields)
    pde1, pde2, pde3, pde4, pde5, pde6, pde7 = pdes
    pre_step_field_update!(pde1, fields)
    pre_step_field_update!(pde2, fields)
    pre_step_field_update!(pde3, fields)
    pre_step_field_update!(pde4, fields)
    pre_step_field_update!(pde5, fields)
    pre_step_field_update!(pde6, fields)
    pre_step_field_update!(pde7, fields)
    return
end

post_step_field_update!(::AbstractPDE, fields) = return

function post_step_update!(pdes, fieldhandler::AbstractFieldHandler, ::BCHandler)
    post_step_field_update!(pdes, fieldhandler.fields)
end

function post_step_field_update!(pdes::TupleOfLength{1, AbstractPDE}, fields)
    pde1, = pdes
    post_step_field_update!(pde1, fields)
    return
end

function post_step_field_update!(pdes::TupleOfLength{2, AbstractPDE}, fields)
    pde1, pde2 = pdes
    post_step_field_update!(pde1, fields)
    post_step_field_update!(pde2, fields)
    return
end

function post_step_field_update!(pdes::TupleOfLength{3, AbstractPDE}, fields)
    pde1, pde2, pde3 = pdes
    post_step_field_update!(pde1, fields)
    post_step_field_update!(pde2, fields)
    post_step_field_update!(pde3, fields)
    return
end

function post_step_field_update!(pdes::TupleOfLength{4, AbstractPDE}, fields)
    pde1, pde2, pde3, pde4 = pdes
    post_step_field_update!(pde1, fields)
    post_step_field_update!(pde2, fields)
    post_step_field_update!(pde3, fields)
    post_step_field_update!(pde4, fields)
    return
end

function post_step_field_update!(pdes::TupleOfLength{5, AbstractPDE}, fields)
    pde1, pde2, pde3, pde4, pde5 = pdes
    post_step_field_update!(pde1, fields)
    post_step_field_update!(pde2, fields)
    post_step_field_update!(pde3, fields)
    post_step_field_update!(pde4, fields)
    post_step_field_update!(pde5, fields)
    return
end

function post_step_field_update!(pdes::TupleOfLength{6, AbstractPDE}, fields)
    pde1, pde2, pde3, pde4, pde5, pde6 = pdes
    post_step_field_update!(pde1, fields)
    post_step_field_update!(pde2, fields)
    post_step_field_update!(pde3, fields)
    post_step_field_update!(pde4, fields)
    post_step_field_update!(pde5, fields)
    post_step_field_update!(pde6, fields)
    return
end

function post_step_field_update!(pdes::TupleOfLength{7, AbstractPDE}, fields)
    pde1, pde2, pde3, pde4, pde5, pde6, pde7 = pdes
    post_step_field_update!(pde1, fields)
    post_step_field_update!(pde2, fields)
    post_step_field_update!(pde3, fields)
    post_step_field_update!(pde4, fields)
    post_step_field_update!(pde5, fields)
    post_step_field_update!(pde6, fields)
    post_step_field_update!(pde7, fields)
    return
end
