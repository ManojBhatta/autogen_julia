export AbstractSolver
export AbstractLinearSolver
export LinearSolver
export solve!
export getsolverstats

import LinearSolve: UMFPACKFactorization
import LinearSolve: KrylovJL_BICGSTAB, KrylovJL_CG, KrylovJL_CRAIGMR, KrylovJL_GMRES, KrylovJL_LSMR, KrylovJL_MINRES
import .LinearSolvePardisoExt: MKLPardisoFactorize, MKLPardisoIterate, PardisoJL

export UMFPACKFactorization
export KrylovJL_BICGSTAB, KrylovJL_CG, KrylovJL_CRAIGMR, KrylovJL_GMRES, KrylovJL_LSMR, KrylovJL_MINRES
export MKLPardisoFactorize, MKLPardisoIterate, PardisoJL

abstract type AbstractSolver end
abstract type AbstractLinearSolver <: AbstractSolver end

function solve!(solver::AbstractLinearSolver, soe, pdes, fieldhandler, bchandler, quadrature=defaultquadrature(fieldhandler); material=UndefinedMaterial())
    applydirichlet!(bchandler, fieldhandler)
    # TODO: Only assemble Jacobian if necessary
    # TODO: Pass quadrature from where?
    assemble!(soe, pdes, fieldhandler, bchandler, material, quadrature)
    solve!(solver, soe)
    update_solution_fields!(pdes, soe)
    write!(solver.writers, pdes, fieldhandler)
    return soe
end

###########################
# LinearSolver
###########################

mutable struct LinearSolver{C, W} <: AbstractLinearSolver
    cache::C
    update_A::Bool
    writers::W
end

# ------------------ LinearSolve.jl ------------------

"""
"""
function LinearSolver(A, b, alg::Union{LinearSolve.SciMLLinearSolveAlgorithm, Nothing}, u0=zero(b);
        update_A=true, writers=nothing, kwargs...)
    condition = LinearSolve.OperatorCondition.IllConditioned
    assumptions = LinearSolve.OperatorAssumptions(true; condition)
    prob = LinearSolve.LinearProblem(A, b; u0, assumptions)
    cache = LinearSolve.init(prob, alg; kwargs...)

    _writers = isnothing(writers) ? nothing : tuple(writers...)

    LinearSolver{typeof(cache), typeof(_writers)}(cache, update_A, _writers)
end

function LinearSolver(soe::SystemOfEquations, alg=LinearSolve.UMFPACKFactorization(); kwargs...)
    LinearSolver(soe.jacobian, soe.residual, alg, soe.Δu; kwargs...)
end

"""
"""
function solve!(solver::LinearSolver{<:LinearSolve.LinearCache}, soe::SystemOfEquations)
    if solver.update_A || solver.cache.isfresh
        solver.cache.A = soe.jacobian
    end
    solver.cache.b = soe.residual
    sol = LinearSolve.solve!(solver.cache)
    return sol
end

function getsolverstats(sol::LinearSolve.SciMLSolution)
    # TODO? Use solver specific stats from sol.cache.cacheval?
    (
        niter = sol.iters,
        # time = stats.timer,
        # ...
    )
end

###########################
#
###########################