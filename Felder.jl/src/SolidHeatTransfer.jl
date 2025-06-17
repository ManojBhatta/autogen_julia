export SolidHeatTransferPDE
export HeatConvectionBC

############################
# PDE Heat Equation
############################

"""
    SolidHeatTransferPDE(field, k, cp, rho, qv=0.0)

Linear heat equation with thermal conductivity `k`, specific heat `cp`, density `rho`
and heat source `qv`.

`k`, `cp`, `rho` and `qv` can be functions of `x` and `t`.

`field` is a `ScalarField` representing the temperature field, which must be
constructed with `cachelength=2` to store the previous time steps for the BDF2 scheme.

At the moment, the time cache must be initialized by hand with something like

    field.t_cache .= [-0.01, 0.0]
"""
@kwdef struct SolidHeatTransferPDE{F, T1, T2, T3, T4} <: AbstractPDE{1, Float64}
    field::F
    k::T1
    cp::T2
    rho::T3
    qv::T4                = 0.0
end

function Felder.weakformvars(pde::SolidHeatTransferPDE, fields, proxy, q)
    T = interpolate(pde.field, q)
    T_old_1 = interpolate_cache(pde.field, 1, q)
    T_old_2 = interpolate_cache(pde.field, 2, q)
    x = proxy.x[q]
    t = pde.field.t[]
    Δt0 = t - pde.field.t_cache[1]

    if Δt0 < Inf
        Δt1 = pde.field.t_cache[1] - pde.field.t_cache[2]
        w = Δt0 / Δt1
        bdf_a = (1 + 2*w) / (1 + w)
        bdf_b = (1 + w)^2 / (1 + w)
        bdf_c = w^2 / (1 + w)
        dTdt = (bdf_a * T - bdf_b * T_old_1 + bdf_c * T_old_2) / Δt0 # BDF2 with variable step
    else
        dTdt = 0.0
        bdf_a = 1.0
    end

    return (
        T = T,
        ∇T = interpolate_grad(pde.field, q),
        k = evalargs(pde.k, x, t),
        cp = evalargs(pde.cp, x, t),
        ρ = evalargs(pde.rho, x, t),
        qv = evalargs(pde.qv, x, t),
        dTdt = dTdt,
        bdf_a = bdf_a,
        Δt = Δt0,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function Felder.residual!(r, ::SolidHeatTransferPDE, vars, proxy, i)
    @unpack T, ∇T, k, cp, ρ, qv, dTdt, bdf_a, Δt = vars
    @unpack N, ∇N = vars
    r[1] = k * ∇T ⋅ ∇N[i] + ρ * cp * dTdt * N[i] - qv * N[i]
end

@inline function Felder.jacobian!(J, ::SolidHeatTransferPDE, vars, proxy, i, j)
    @unpack T, ∇T, k, cp, ρ, qv, dTdt, bdf_a, Δt = vars
    @unpack N, ∇N = vars
    J[1, 1] = k * ∇N[j] ⋅ ∇N[i] + ρ * cp * bdf_a * N[j] / Δt * N[i]
end

function post_step_field_update!(pde::SolidHeatTransferPDE, fields)
    pde.field.u_cache[2] .= pde.field.u_cache[1]
    pde.field.u_cache[1] .= pde.field.u
    pde.field.t_cache[2] = pde.field.t_cache[1]
    pde.field.t_cache[1] = pde.field.t[]
end

########################
# %% Heat Convection BC
########################

struct HeatConvectionBC{F1, F2, F3} <: AbstractIntegratedBC
    field::F3
    boundaries::Vector{Symbol}
    h::F1 # h(x, t)
    T0::F2 # T0(x, t)
end

function HeatConvectionBC(field, boundarytags, h, T0)
    _btags = boundarytags isa Symbol ? [boundarytags] : boundarytags

    HeatConvectionBC{typeof(a), typeof(b), typeof(field)}(field, _btags, h, T0)
end

function Base.show(io::Core.IO, ::MIME"text/plain", bc::HeatConvectionBC)
    println(io, "HeatConvectionBC for:")
    println(io, "  Field:            $(typeof(bc.field).name.wrapper){$(eltype(bc.field)), ...}")
    println(io, "  Boundaries:       $(bc.boundaries)")
    println(io, "  h:                $(bc.h)")
    print(io,   "  T0:                $(bc.T0)")
end

function bcvars(bc::HeatConvectionBC, fields, proxy, q)
    x = proxy.x[q]
    t = bc.field.t[]
    T = interpolate(bc.field, q)
    T0 = evalargs(bc.T0, x, t)
    return (
        ΔT = T - T0,
        h = evalargs(bc.h, x, t),
        x = x,
        t = t,
        N = proxy.N[q],
    )
end

@inline function residual!(r, ::HeatConvectionBC, vars, proxy, i)
    @unpack h, ΔT, N = vars
    r[1] = h * ΔT * N[i]
end

@inline function jacobian!(J, ::HeatConvectionBC, vars, proxy, i, j)
    @unpack h, ΔT, N = vars
    J[1, 1] = h * N[j] * N[i]
end
