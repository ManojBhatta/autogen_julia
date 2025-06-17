export AbstractMaterial
export UndefinedMaterial
export AbstractMaterialHandler
export MaterialHandler

abstract type AbstractMaterial end

struct UndefinedMaterial <: AbstractMaterial end
