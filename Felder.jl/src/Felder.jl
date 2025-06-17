module Felder

using StaticArrays
using LinearAlgebra
using Printf
using SparseArrays
using OrderedCollections
using WriteVTK
using Parameters: @unpack, @with_kw
using NearestNeighbors
using Base.Threads: @threads, nthreads, threadid
using Base: RefValue, @kwdef, @propagate_inbounds
# using MKLSparse # Does not work together with Pardiso

import Krylov
import LinearSolve
include("LinearSolvePardisoExt.jl")
import .LinearSolvePardisoExt
import Preconditioners
using Preconditioners: CholeskyPreconditioner
import IncompleteLU

export @unpack

# Alias for Tuple or NamedTuple of length N and type T
const TupleOfLength{N, T} = Union{NamedTuple{<:Any, <:NTuple{N, T}}, NTuple{N, T}}

include("Misc.jl")
include("Shapes.jl")
include("Mesh.jl")
include("MeshGeneration.jl")
include("UnvMesh.jl")
include("ShapeFunctions.jl")
include("Lagrange.jl")
include("Quadrature.jl")
include("Material.jl")
include("DofMap.jl")
include("Proxy.jl")
include("Fields.jl")
include("PDE.jl")
include("BCs.jl")
include("SystemOfEquations.jl")
include("Assembly.jl")
include("Probes.jl")
include("Distance.jl")
include("FieldInterpolation.jl")
include("FieldIntegration.jl")
include("LinearSolver.jl")
include("NonlinearSolver.jl")
include("TransientSolver.jl")
include("VTKOutput.jl")

# Physics Implementations
include("SolidHeatTransfer.jl")

# For overloading in ../ext/FelderMakieExt.jl
export meshlineplot, meshlineplot!
export meshplot, meshplot!
export fieldlineplot, fieldlineplot!
export fieldplot, fieldplot!
export fieldsurface, fieldsurface!
export to_gb_mesh

function meshlineplot end
function meshlineplot! end
function meshplot end
function meshplot! end
function fieldlineplot end
function fieldlineplot! end
function fieldplot end
function fieldplot! end
function fieldsurface end
function fieldsurface! end
function to_gb_mesh end

end # module
