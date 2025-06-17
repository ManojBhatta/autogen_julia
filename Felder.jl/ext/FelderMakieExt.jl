module FelderMakieExt

using Felder
using StaticArrays
using LinearAlgebra
import Makie
import Makie: GeometryBasics

import Felder: meshlineplot, meshlineplot!
import Felder: meshplot, meshplot!
import Felder: fieldlineplot, fieldlineplot!
import Felder: fieldplot, fieldplot!
import Felder: fieldsurface, fieldsurface!
import Felder: to_gb_mesh

include("FelderMakieRecipes.jl")

end # module