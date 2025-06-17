using Felder
using LinearAlgebra
using StaticArrays
import Zygote
import Felder: weakformvars, residual!, jacobian!, post_step_field_update!, evalargs, bcvars
import Felder: nonlinear_field_update!
using GLMakie: GLMakie, Makie, Figure, Axis, Axis3, plot, plot!, scatter!
using GLMakie: Observable, Label, Menu, Toggle, Button, Textbox, GridLayout
using Printf

########################
# %% Physics Definition
########################

@kwdef struct CoefficientFormPDE1D{F, T1, T2, T3, T4, T5, T6} <: AbstractPDE{1, Float64}
    field::F
    delta::T1
    nu::T2
    gamma::T3
    beta::T4
    alpha::T5
    g::T6
    sleeptime::Float64 = 0.0
end

function weakformvars(pde::CoefficientFormPDE1D, fields, proxy, q)
    u = interpolate(pde.field, q)
    u_old = interpolate_cache(pde.field, 1, q)
    ∇u = interpolate_grad(pde.field, q)[1]
    x = proxy.x[q][1]
    t = pde.field.t[]
    Δt = t - pde.field.t_cache[1]

    dudt = (u - u_old) / Δt # BDF1
    ddudt_du = 1 / Δt

    delta, (ddelta_du, ddelta_d∇u, ) = Zygote.withgradient(pde.delta, u, ∇u, x, t)
    nu, (dnu_du, dnu_d∇u, ) = Zygote.withgradient(pde.nu, u, ∇u, x, t)
    gamma, (dgamma_du, dgamma_d∇u, ) = Zygote.withgradient(pde.gamma, u, ∇u, x, t)
    beta, (dbeta_du, dbeta_d∇u, ) = Zygote.withgradient(pde.beta, u, ∇u, x, t)
    alpha, (dalpha_du, dalpha_d∇u, ) = Zygote.withgradient(pde.alpha, u, ∇u, x, t)
    g, (dg_du, dg_d∇u, ) = Zygote.withgradient(pde.g, u, ∇u, x, t)

    nu_nl = nu
    isnothing(dnu_d∇u)    || (nu_nl += dnu_d∇u * ∇u)
    isnothing(dgamma_d∇u) || (nu_nl += dgamma_d∇u * u)

    gamma_nl = gamma
    isnothing(dnu_du)    || (gamma_nl += dnu_du * ∇u)
    isnothing(dgamma_du) || (gamma_nl += dgamma_du * u)

    beta_nl = beta
    isnothing(dbeta_d∇u)  || (beta_nl += dbeta_d∇u * ∇u)
    isnothing(dalpha_d∇u) || (beta_nl += dalpha_d∇u * u)
    isnothing(dg_d∇u)     || (beta_nl -= dg_d∇u)
    isnothing(ddelta_d∇u) || (beta_nl += ddelta_d∇u * dudt)

    alpha_nl = alpha
    isnothing(dbeta_du)  || (alpha_nl += dbeta_du * ∇u)
    isnothing(dalpha_du) || (alpha_nl += dalpha_du * u)
    isnothing(dg_du)     || (alpha_nl -= dg_du)
    alpha_nl += delta * ddudt_du
    isnothing(ddelta_du) || (alpha_nl += ddelta_du * dudt)

    return (
        u = u,
        ∇u = ∇u,
        dudt = dudt,
        delta = delta,
        nu = nu,
        gamma = gamma,
        beta = beta,
        alpha = alpha,
        g = g,
        nu_nl = nu_nl,
        gamma_nl = gamma_nl,
        beta_nl = beta_nl,
        alpha_nl = alpha_nl,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(r, ::CoefficientFormPDE1D, vars, proxy, i)
    @unpack u, ∇u, dudt, N, ∇N = vars
    @unpack delta, nu, gamma, beta, alpha, g = vars

    r[1] = (nu * ∇u + gamma * u) * ∇N[i][1] + (beta * ∇u + alpha * u - g + delta * dudt) * N[i]
end

@inline function jacobian!(J, ::CoefficientFormPDE1D, vars, proxy, i, j)
    @unpack u, ∇u, N, ∇N = vars
    @unpack nu_nl, gamma_nl, beta_nl, alpha_nl= vars

    J[1, 1] = (nu_nl * ∇N[j][1] + gamma_nl * N[j]) * ∇N[i][1] + (beta_nl ⋅ ∇N[j][1] + alpha_nl * N[j]) * N[i]
end

function nonlinear_field_update!(::AbstractPDE, fields)
    if isabort[]
        throw(InterruptException())
    end
end

function post_step_field_update!(pde::CoefficientFormPDE1D, fields)
    fields.u.u_cache[1] .= fields.u.u
    fields.u.t_cache[1] = fields.u.t[]
    Makie.notify(field_node)
    sleep(pde.sleeptime) # For animation
end

########################
# %% BC Definition
########################

struct GeneralBC{F1, F2, F3, F4} <: AbstractIntegratedBC
    field::F4
    boundaries::Vector{Symbol}
    a::F1 # a(u, ux, x, t)
    b::F2 # b(u, ux, x, t)
    c::F3 # c(u, ux, x, t)
end

function bcvars(bc::GeneralBC, fields, proxy, q)
    u = interpolate(bc.field, q)
    ∇u = interpolate_grad(bc.field, q)[1]
    x = proxy.x[q][1]
    t = bc.field.t[]
    a, (da_du, da_d∇u, ) = Zygote.withgradient(bc.a, u, ∇u, x, t)
    b, (db_du, db_d∇u, ) = Zygote.withgradient(bc.b, u, ∇u, x, t)
    c, (dc_du, dc_d∇u, ) = Zygote.withgradient(bc.c, u, ∇u, x, t)

    c_nl = c
    isnothing(dc_d∇u) || (c_nl += dc_d∇u * ∇u)
    isnothing(db_d∇u) || (c_nl += db_d∇u * u)
    isnothing(da_d∇u) || (c_nl += da_d∇u)

    b_nl = b
    isnothing(dc_du) || (b_nl += dc_du * ∇u)
    isnothing(db_du) || (b_nl += db_du * u)
    isnothing(da_du) || (b_nl += da_du)

    return (
        u = u,
        ∇u = ∇u,
        a = a,
        b = b,
        c = c,
        b_nl = b_nl,
        c_nl = c_nl,
        N = proxy.N[q],
        ∇N = proxy.dNdx[q],
    )
end

@inline function residual!(r, ::GeneralBC, vars, proxy, i)
    @unpack a, b, c, u, ∇u, N, ∇N = vars
    r[1] = (c * ∇u + b * u + a) * N[i]
end

@inline function jacobian!(J, ::GeneralBC, vars, proxy, i, j)
    @unpack b_nl, c_nl, N, ∇N = vars
    J[1, 1] = (c_nl * ∇N[j][1] + b_nl * N[j]) * N[i]
end

########################
# Presets
########################

presets = Dict(
    "Zero defaults" => (
        delta = "0",
        nu = "0",
        gamma = "0",
        beta = "0",
        alpha = "0",
        g = "0",
        bc_left = Nothing,
        c_left = "0",
        b_left = "0",
        a_left = "0",
        bc_right = Nothing,
        c_right = "0",
        b_right = "0",
        a_right = "0",
        u_init = "0",
        dudt_init = "0",
        sf = 2,
        mesh = "range(-1, 1, step=0.1)",
        timesteps = "range(0, 2, length=100)",
        quadorder = "5",
        damp = "1.0",
        sleep = "0.0",
    ),
    "Example Expressions" => (
        delta = "1",
        nu = "x < 5 ? exp(u) : 1",
        gamma = "1e-2 * u^2",
        beta = "u",
        alpha = "-ux",
        g = "t < 5 ? cos(2 * pi * x + t) : 0",
        bc_left = DirichletBC,
        c_left = "0",
        b_left = "0",
        a_left = "t / 5",
        bc_right = GeneralBC,
        c_right = "x < 5 ? exp(u) : 1",
        b_right = "-1e-2 * u^2",
        a_right = "0",
        u_init = "x < 5 ? (5 - x) : (5 * sin(pi * x))",
        dudt_init = "min(0, x)",
        sf = 2,
        mesh = "10 .^ range(0, 1, length=40) .- 1",
        timesteps = "[0, 0.1, 0.2, 0.5, 1, 2]",
    ),
    "Heat Equation" => (
        delta = "1.0",
        nu = "1.0",
        gamma = "0",
        beta = "0",
        alpha = "0",
        g = "0",
        bc_left = Nothing,
        c_left = "0",
        b_left = "0",
        a_left = "0",
        bc_right = Nothing,
        c_right = "0",
        b_right = "0",
        a_right = "0",
        u_init = "abs(x) > 0.5 ? 0 : 1",
        dudt_init = "0",
        mesh = "range(-1, 1, step=0.1)",
        timesteps = "range(0, 1, length=100)",
    ),
    "Heat Equation 2" => (
        delta = "1.0",
        nu = "1.0",
        gamma = "0",
        beta = "0",
        alpha = "0",
        g = "1.0",
        bc_left = DirichletBC,
        c_left = "0",
        b_left = "0",
        a_left = "0.0",
        bc_right = DirichletBC,
        c_right = "0",
        b_right = "0",
        a_right = "0.0",
        u_init = "0",
        dudt_init = "0",
        mesh = "range(-1, 1, step=0.1)",
        timesteps = "range(0, 2, length=100)",
    ),
    "Fisher's Equation" => (
        delta = "1.0",
        nu = "1.0",
        gamma = "0",
        beta = "0",
        alpha = "u - 1",
        g = "0",
        bc_left = Nothing,
        c_left = "0",
        b_left = "0",
        a_left = "0",
        bc_right = Nothing,
        c_right = "0",
        b_right = "0",
        a_right = "0",
        u_init = "abs(x) < 5 ? 0.1 : 0",
        dudt_init = "0",
        mesh = "range(-60, 60, length=120)",
        timesteps = "range(0, 30, length=200)",
    ),
    "Approximator" => (
        delta = "0",
        nu = "0",
        gamma = "0",
        beta = "0",
        alpha = "1.0",
        g = "sin(2 * pi * x * t)",
        bc_left = Nothing,
        c_left = "0",
        b_left = "0",
        a_left = "0",
        bc_right = Nothing,
        c_right = "0",
        b_right = "0",
        a_right = "0",
        u_init = "0",
        dudt_init = "0",
        mesh = "range(-1, 1, step=0.1)",
        timesteps = "range(0, 3, length=200)",
        damp = "1.0",
    ),
    "Convection-Diffusion" => (
        delta = "1.0",
        nu = "0.05",
        gamma = "-1.0",
        beta = "0",
        alpha = "0",
        g = "abs(x) < 0.2 ? sign(sinpi(t)) : 0",
        bc_left = DirichletBC,
        c_left = "0",
        b_left = "0",
        a_left = "0.0",
        bc_right = GeneralBC,
        c_right = "0.05",
        b_right = "1.0",
        a_right = "0.0",
        u_init = "0",
        dudt_init = "0",
        mesh = "range(-1, 5, step=0.02)",
        timesteps = "range(0, 10, length=220)",
        damp = "1.0",
    ),
    "Convection-Diffusion Test" => (
        delta = "1.0",
        nu = "0.05",
        gamma = "-1.0",
        beta = "0",
        alpha = "0",
        g = "0",
        bc_left = DirichletBC,
        c_left = "0",
        b_left = "0",
        a_left = "1.0",
        bc_right = DirichletBC,
        c_right = "0",
        b_right = "0",
        a_right = "0.0",
        u_init = "0",
        dudt_init = "0",
        mesh = "range(-1, 1, step=0.05)",
        timesteps = "range(0, 10, length=100)",
        damp = "1.0",
    ),
    "Burgers' Equation" => (
        delta = "1.0",
        nu = "1.0",
        gamma = "0",
        beta = "u",
        alpha = "0",
        g = "0",
        bc_left = Nothing,
        c_left = "0",
        b_left = "0",
        a_left = "0",
        bc_right = Nothing,
        c_right = "0",
        b_right = "0",
        a_right = "0",
        u_init = "exp(-x^2)",
        dudt_init = "0",
        mesh = "range(-3, 3, step=0.01)",
        timesteps = "range(0, 3, length=500)",
    ),
    "Burgers' Equation (inviscid)" => (
        delta = "1.0",
        nu = "0",
        gamma = "0",
        beta = "u",
        alpha = "0",
        g = "0",
        bc_left = Nothing,
        c_left = "0",
        b_left = "0",
        a_left = "0",
        bc_right = Nothing,
        c_right = "0",
        b_right = "0",
        a_right = "0",
        u_init = "exp(-x^2)",
        dudt_init = "0",
        mesh = "range(-3, 3, step=0.01)",
        timesteps = "range(0, 3, length=500)",
    ),
    "Fujita-Storm Equation" => (
        delta = "3.0",
        nu = "u^-2",
        gamma = "0",
        beta = "0",
        alpha = "0",
        g = "0",
        bc_left = Nothing,
        c_left = "0",
        b_left = "0",
        a_left = "0",
        bc_right = Nothing,
        c_right = "0",
        b_right = "0",
        a_right = "0",
        u_init = "exp(-x^2)",
        dudt_init = "0",
        mesh = "range(-3, 3, step=0.1)",
        timesteps = "range(0, 1, length=500)",
    ),
    "Poisson Equation" => (
        delta = "0",
        nu = "1",
        gamma = "0",
        beta = "0",
        alpha = "0",
        g = "abs(x) < t ? 1 : 0",
        bc_left = DirichletBC,
        c_left = "0",
        b_left = "0",
        a_left = "0",
        bc_right = Nothing,
        c_right = "0",
        b_right = "0",
        a_right = "0",
        u_init = "0",
        dudt_init = "0",
        mesh = "range(-1, 1, step=0.05)",
        timesteps = "range(0, 0.5, length=100)",
    ),
)

########################
# GUI Elements
########################

function quad_validator(str)
    try
        i = parse(Int64, str)
        if 0 < i <= 40
            return true
        end
    catch
        return false
    end
    return false
end

isbusy = Ref(false)
isabort = Ref(false)

GLMakie.activate!(;
    float=false, # Always float on top
    focus_on_show=true, # Focus the window when newly opened
    decorated=true # Show window decorations
)

Makie.set_theme!(Makie.theme_dark())
Makie.update_theme!(backgroundcolor="#0D1117") # GitHub Dark Theme Background

fig = Makie.Figure(size=(1280, 760), fontsize=14)

grid_A = fig[1, 1, Makie.TopLeft()] = GridLayout(width=350)
grid_B = fig[1, 2, Makie.TopLeft()] = GridLayout()
grid_B1 = grid_B[1, 1] = GridLayout()
grid_B2 = grid_B[2, 1] = GridLayout()

grid_settings = grid_A

grid_pde = grid_B1[1, :] = GridLayout()
grid_coeffs = grid_B1[2, :] = GridLayout()
grid_bc = grid_B1[3, :] = GridLayout()

grid_bc_left = grid_bc[:, 1] = GridLayout()
grid_bc_right = grid_bc[:, 2] = GridLayout()

# --------------------------------

Label(grid_pde[1, 1], Makie.@L_str("\\delta \\dot{u} - \\mathbf{\\nabla} \\cdot (\\nu \\mathbf{\\nabla} u + \\mathbf{\\gamma} u) + \\mathbf{\\beta} \\cdot \\mathbf{\\nabla} u + \\alpha u = g"),
    font=:regular, tellheight=true, fontsize=24)

# --------------------------------

Label(grid_coeffs[1, 1], Makie.@L_str("\\delta(u, u_x, x, t)"), fontsize=20)
tb_delta = Textbox(grid_coeffs[2, 1], stored_string="1.0")

Label(grid_coeffs[1, 2], Makie.@L_str("\\nu(u, u_x, x, t)"), fontsize=20)
tb_nu = Textbox(grid_coeffs[2, 2], stored_string="0.1")

Label(grid_coeffs[1, 3], Makie.@L_str("\\gamma(u, u_x, x, t)"), fontsize=20)
tb_gamma = Textbox(grid_coeffs[2, 3], stored_string="0.0")

Label(grid_coeffs[1, 4], Makie.@L_str("\\beta(u, u_x, x, t)"), fontsize=20)
tb_beta = Textbox(grid_coeffs[2, 4], stored_string="0.0")

Label(grid_coeffs[1, 5], Makie.@L_str("\\alpha(u, u_x, x, t)"), fontsize=20)
tb_alpha = Textbox(grid_coeffs[2, 5], stored_string="0.0")

Label(grid_coeffs[1, 6], Makie.@L_str("g(u, u_x, x, t)"), fontsize=20)
tb_g = Textbox(grid_coeffs[2, 6], stored_string="0.0")


# --------------------------------

bc_options = [
    ("None", Nothing),
    ("Dirichet u = a(x, t)", DirichletBC),
    ("Contribution -(ν∇u + γu)|ₙ = -(c∇u + bu + a)|ₙ", GeneralBC),
]

Label(grid_bc_left[1, 1][1, 1], "Left boundary:", halign=:left)
menu_bc_left = Menu(grid_bc_left[1, 1][1, 2], options=bc_options, default=2, width=300, halign=:left)
label_left_c = Label(grid_bc_left[2, 1][1, 1], Makie.@L_str("c(u, u_x, x, t)"), fontsize=20, halign=:center)
tb_bc_left_c = Textbox(grid_bc_left[2, 1][2, 1], stored_string="0.0", halign=:center)
label_left_b = Label(grid_bc_left[2, 1][1, 2], Makie.@L_str("b(u, u_x, x, t)"), fontsize=20, halign=:center)
tb_bc_left_b = Textbox(grid_bc_left[2, 1][2, 2], stored_string="0.0", halign=:center)
label_left_a = Label(grid_bc_left[2, 1][1, 3], Makie.@L_str("a(u, u_x, x, t)"), fontsize=20, halign=:center)
tb_bc_left_a = Textbox(grid_bc_left[2, 1][2, 3], stored_string="cos(2 * pi * t)", halign=:center)

Label(grid_bc_right[1, 1][1, 1], "Right boundary:", halign=:left)
menu_bc_right = Menu(grid_bc_right[1, 1][1, 2], options=bc_options, default=2, width=300, halign=:left)
label_right_c = Label(grid_bc_right[2, 1][1, 1], Makie.@L_str("c(u, u_x, x, t)"), fontsize=20, halign=:center)
tb_bc_right_c = Textbox(grid_bc_right[2, 1][2, 1], stored_string="0.0", halign=:center)
label_right_b = Label(grid_bc_right[2, 1][1, 2], Makie.@L_str("b(u, u_x, x, t)"), fontsize=20, halign=:center)
tb_bc_right_b = Textbox(grid_bc_right[2, 1][2, 2], stored_string="0.0", halign=:center)
label_right_a = Label(grid_bc_right[2, 1][1, 3], Makie.@L_str("a(u, u_x, x, t)"), fontsize=20, halign=:center)
tb_bc_right_a = Textbox(grid_bc_right[2, 1][2, 3], stored_string="cos(2 * pi * t)", halign=:center)

# --------------------------------

Label(grid_settings[1, 1:2],
    """Solving the weak form of the following PDE in 1D:
    Time integrator: BDF1
    Linear solver: UMFPACKFactorization""", justification=:left,  halign=:left, font=:regular)


sf_options = [
    ("Lagrange 1", Lagrange(1)),
    ("Lagrange 2", Lagrange(2)),
    ("Lagrange 3", Lagrange(3)),
]

preset_options = sort(collect(keys(presets)))

Label(grid_settings[2, 1], "Preset:", halign=:left)
menu_preset = Menu(grid_settings[2, 2], options=preset_options, default=3, width=200, halign=:left)

Label(grid_settings[3, 1], Makie.@L_str("u_\\text{init}(x, t) ="), fontsize=20, halign=:left)
tb_u_init = Textbox(grid_settings[3, 2], stored_string="0.0", halign=:left)

Label(grid_settings[4, 1], Makie.@L_str("\\dot{u}_\\text{init}(x, t) ="), fontsize=20, halign=:left)
tb_dudt_init = Textbox(grid_settings[4, 2], stored_string="0.0", halign=:left)

Label(grid_settings[5, 1], "Discretization:", halign=:left)
menu_sf = Menu(grid_settings[5, 2], options=sf_options, default=2,
    width=100, halign=:left
)

Label(grid_settings[6, 1], "Mesh x:", halign=:left)
tb_mesh = Textbox(grid_settings[6, 2], stored_string="range(-1, 1, step=0.1)", halign=:left)
Label(grid_settings[7, 2][1, 1], "Plot nodes", halign=:left)
tg_nodes = Toggle(grid_settings[7, 2][1, 2], active=false, halign=:left)

Label(grid_settings[8, 1], "Time steps t:", halign=:left)
tb_timesteps = Textbox(grid_settings[8, 2], stored_string="range(0, 2, length=100)", halign=:left)

Label(grid_settings[9, 1], "Quadrature order:", halign=:left)
tb_quadorder = Textbox(grid_settings[9, 2], stored_string="5", validator=quad_validator, halign=:left)

Label(grid_settings[10, 1], "Newton damping:", halign=:left)
tb_damp = Textbox(grid_settings[10, 2], stored_string="1.0", validator=Float64, halign=:left)

Label(grid_settings[11, 1], "Frame sleep [s]:", halign=:left)
tb_sleep = Textbox(grid_settings[11, 2], stored_string="0.0", validator=Float64, halign=:left)

Label(grid_settings[12, 1], "Axis limits:", halign=:left)
Label(grid_settings[12, 2][1, 1], "Auto", halign=:left)
tg_axlimits = Toggle(grid_settings[12, 2][1, 2], active=true, halign=:left)
tb_axlimit_min = Textbox(grid_settings[13, 2][1, 1], stored_string="-1.0", validator=Float64, halign=:left)
Label(grid_settings[13, 2][1, 2], "to", halign=:left)
tb_axlimit_max = Textbox(grid_settings[13, 2][1, 3], stored_string="1.0", validator=Float64, halign=:left)

Label(grid_settings[14, 1], "Time history length:", halign=:left)
tb_history_length = Textbox(grid_settings[14, 2][1, 1], stored_string="25", validator=Int, halign=:left)
Label(grid_settings[14, 2][1, 2], "Plot", halign=:left)
tg_history = Toggle(grid_settings[14, 2][1, 3], active=false, halign=:left)

# --------------------------------

button_calculate = Button(grid_settings[15, 1:2], label="Calculate", labelcolor=:black, width=120)

# --------------------------------

axtitle = Observable("t = ")
ax1 = Axis(grid_B2[1, 1], title=axtitle, xlabel="x", ylabel="u", alignmode=Makie.Outside(), xautolimitmargin=(0, 0))

########################
# Interaction
########################

evalparse(tb::Textbox) = eval(Meta.parse(tb.stored_string[]))

Makie.on(button_calculate.clicks) do clicks
    if isbusy[]
        isabort[] = true
        return
    end

    fields, model... = getmodel()

    global field_node = Observable(fields.u)
    Makie.on(x -> (axtitle[] = string("t = ",  @sprintf("%g", x.t[]))), field_node)

    global fieldhistory = Observable([deepcopy(fields.u)])

    history_stride = max(1, length(model[1].timesteps) ÷ parse(Int, tb_history_length.stored_string[]))
    counter = 0
    Makie.on(field_node) do field
        counter += 1
        if counter % history_stride == 0
            push!(fieldhistory[], deepcopy(field))
        end
    end

    empty!(ax1)

    tg_nodes.active.val = false
    pl = plot!(ax1, field_node, color=Makie.wong_colors()[1], pointsperelement=max(2, 500÷nelements(fields.u.mesh)))

    if tg_axlimits.active[]
        ymin = minimum(fields.u.u)
        ax1.yaxis.attributes.limits.val = (ymin, ymin + 10 * eps())
        Makie.on(x -> grow_axis_ylimits!(ax1), field_node)
    end

    Makie.@async begin
        isbusy[] = true
        button_calculate.label[] = "Abort"
        try
            solve!(model...)
            display(model[2])
        catch err
            showerror(stdout, err, catch_backtrace())
            println("Abort.")
        end
        update_limit_textboxes(ax1)
        isbusy[] = false
        isabort[] = false
        button_calculate.label[] = "Calculate"

        if tg_history.active[]
            plot_history()
        end
    end
end

Makie.on(tg_history.active, priority=0) do active
    if active
        plot_history()
    else
        for pl in pl_history
            delete!(ax1, pl)
        end
    end
end

Makie.on(tb_delta.stored_string) do stored_string
    eval(Meta.parse("delta(u, ux, x, t) = " * stored_string))
    println("Update delta(u, ux, x, t) = $stored_string")
end

Makie.on(tb_nu.stored_string) do stored_string
    eval(Meta.parse("nu(u, ux, x, t) = " * stored_string))
    println("Update nu(u, ux, x, t) = $stored_string")
end

Makie.on(tb_gamma.stored_string) do stored_string
    eval(Meta.parse("gamma(u, ux, x, t) = " * stored_string))
    println("Update gamma(u, ux, x, t) = $stored_string")
end

Makie.on(tb_beta.stored_string) do stored_string
    eval(Meta.parse("beta(u, ux, x, t) = " * stored_string))
    println("Update beta(u, ux, x, t) = $stored_string")
end

Makie.on(tb_alpha.stored_string) do stored_string
    eval(Meta.parse("alpha(u, ux, x, t) = " * stored_string))
    println("Update alpha(u, ux, x, t) = $stored_string")
end

Makie.on(tb_g.stored_string) do stored_string
    eval(Meta.parse("g(u, ux, x, t) = " * stored_string))
    println("Update g(u, ux, x, t) = $stored_string")
end

Makie.on(tb_u_init.stored_string, priority=1) do stored_string
    eval(Meta.parse("u_init(x, t) = " * stored_string))
    println("Update u_init(x, t) = $stored_string")
end

Makie.on(tb_u_init.stored_string, priority=0) do stored_string
    plot_u_init()
end

Makie.on(tb_dudt_init.stored_string) do stored_string
    eval(Meta.parse("dudt_init(x, t) = " * stored_string))
    println("Update dudt_init(x, t) = $stored_string")
end

Makie.on(tb_bc_left_a.stored_string) do stored_string
    eval(Meta.parse("a_left(u, ux, x, t) = " * stored_string))
    println("Update a_left(u, ux, x, t) = $stored_string")
end

Makie.on(tb_bc_left_b.stored_string) do stored_string
    eval(Meta.parse("b_left(u, ux, x, t) = " * stored_string))
    println("Update b_left(u, ux, x, t) = $stored_string")
end

Makie.on(tb_bc_left_c.stored_string) do stored_string
    eval(Meta.parse("c_left(u, ux, x, t) = " * stored_string))
    println("Update c_left(u, ux, x, t) = $stored_string")
end

Makie.on(tb_bc_right_a.stored_string) do stored_string
    eval(Meta.parse("a_right(u, ux, x, t) = " * stored_string))
    println("Update a_right(u, ux, x, t) = $stored_string")
end

Makie.on(tb_bc_right_b.stored_string) do stored_string
    eval(Meta.parse("b_right(u, ux, x, t) = " * stored_string))
    println("Update b_right(u, ux, x, t) = $stored_string")
end

Makie.on(tb_bc_right_c.stored_string) do stored_string
    eval(Meta.parse("c_right(u, ux, x, t) = " * stored_string))
    println("Update c_right(u, ux, x, t) = $stored_string")
end

Makie.on(menu_bc_left.selection) do selection
    println("Select left BC: $selection")
    if selection == Nothing
        label_left_c.text[] = "       -       "
        label_left_b.text[] = "       -       "
        label_left_a.text[] = "       -       "
    elseif selection == DirichletBC
        label_left_c.text[] = "       -       "
        label_left_b.text[] = "       -       "
        label_left_a.text[] = Makie.@L_str("a(x, t)")
    else
        label_left_c.text[] = Makie.@L_str("c(u, u_x, x, t)")
        label_left_b.text[] = Makie.@L_str("b(u, u_x, x, t)")
        label_left_a.text[] = Makie.@L_str("a(u, u_x, x, t)")
    end
end

Makie.on(menu_bc_right.selection) do selection
    println("Select right BC: $selection")
    if selection == Nothing
        label_right_c.text[] = "       -       "
        label_right_b.text[] = "       -       "
        label_right_a.text[] = "       -       "
    elseif selection == DirichletBC
        label_right_c.text[] = "       -       "
        label_right_b.text[] = "       -       "
        label_right_a.text[] = Makie.@L_str("a(x, t)")
    else
        label_right_c.text[] = Makie.@L_str("c(u, u_x, x, t)")
        label_right_b.text[] = Makie.@L_str("b(u, u_x, x, t)")
        label_right_a.text[] = Makie.@L_str("a(u, u_x, x, t)")
    end
end

Makie.on(menu_sf.selection) do selection
    println("Select discretization: $selection")
    plot_u_init()
end

Makie.on(tb_mesh.stored_string) do stored_string
    println("Set mesh: $stored_string")
    plot_u_init()
end

Makie.on(tg_nodes.active) do active
    println("Set plot nodes: $active")
    if active
        plot_nodes()
    else
        Makie.delete!(ax1, pl_nodes)
    end
end

Makie.on(tb_timesteps.stored_string) do stored_string
    println("Set timesteps: $stored_string")
    plot_u_init()
end

Makie.on(tg_axlimits.active) do active
    if active
        println("Set axis limits: auto")
        Makie.autolimits!(ax1)
        update_limit_textboxes(ax1)
    else
        println("Set axis limits: manual")
        ylim_min = parse(Float64, tb_axlimit_min.stored_string[])
        ylim_max = parse(Float64, tb_axlimit_max.stored_string[])
        Makie.ylims!(ax1, (ylim_min, ylim_max))
    end
end

Makie.on(tb_axlimit_min.stored_string) do stored_string
    tg_axlimits.active[] = false
    println("Set axis limit ymin: $stored_string")
end

Makie.on(tb_axlimit_max.stored_string) do stored_string
    tg_axlimits.active[] = false
    println("Set axis limit ymax: $stored_string")
end

function plot_u_init()
    fields, model... = getmodel()

    global field_node = Observable(fields.u)
    axtitle[] = string("t = ",  @sprintf("%g", fields.u.t[]))

    global fieldhistory = Observable(typeof(fields.u)[])

    Makie.empty!(ax1)
    plot!(ax1, field_node, color=Makie.wong_colors()[1], pointsperelement=max(2, 500÷nelements(fields.u.mesh)))

    if tg_nodes.active[]
        plot_nodes()
    end

    if tg_axlimits.active[]
        Makie.autolimits!(ax1)
        update_limit_textboxes(ax1)
    end
end

function plot_nodes()
    global pl_nodes = plot!(ax1, field_node, color=Makie.wong_colors()[1], pointsperelement=2, linewidth=0, markersize=9)
end

function plot_history()
    if length(fieldhistory[]) == 0
        global pl_history = typeof(field_node[])[]
        return
    end
    colors = Makie.cgrad(:viridis, length(fieldhistory[]), categorical=true)
    global pl_history = [
        plot!(ax1, field, color=color, pointsperelement=max(2, 500÷nelements(field.mesh)))
        for (field, color) in zip(fieldhistory[], colors)
    ]
end

function grow_axis_ylimits!(ax1)
    ylimits = ax1.yaxis.attributes.limits[]
    ymin, ymax = Makie.autolimits(ax1, 2)
    Makie.ylims!(ax1, (min(ylimits[1], ymin), max(ylimits[2], ymax)))
    return ax1
end

function update_limit_textboxes(ax1)
    tb_axlimit_min.displayed_string[] = @sprintf("%.3g", ax1.yaxis.attributes.limits[][1])
    tb_axlimit_max.displayed_string[] = @sprintf("%.3g", ax1.yaxis.attributes.limits[][2])
    tb_axlimit_min.stored_string.val = tb_axlimit_min.displayed_string[]
    tb_axlimit_max.stored_string.val = tb_axlimit_max.displayed_string[]
end

# ----------------------------------

Makie.on(menu_preset.selection) do selection
    println("Select preset: $selection")
    apply_preset(presets[selection])
end

function apply_preset(preset)
    apply_preset!(tb_delta, preset, :delta)
    apply_preset!(tb_nu, preset, :nu)
    apply_preset!(tb_gamma, preset, :gamma)
    apply_preset!(tb_beta, preset, :beta)
    apply_preset!(tb_alpha, preset, :alpha)
    apply_preset!(tb_g, preset, :g)

    apply_bc_preset!(menu_bc_left, preset, :bc_left)
    apply_bc_preset!(menu_bc_right, preset, :bc_right)
    apply_preset!(tb_bc_left_a, preset, :a_left)
    apply_preset!(tb_bc_left_b, preset, :b_left)
    apply_preset!(tb_bc_left_c, preset, :c_left)
    apply_preset!(tb_bc_right_a, preset, :a_right)
    apply_preset!(tb_bc_right_b, preset, :b_right)
    apply_preset!(tb_bc_right_c, preset, :c_right)

    apply_preset!(tb_u_init, preset, :u_init)
    apply_preset!(tb_dudt_init, preset, :dudt_init)
    :sf in keys(preset) && (menu_sf.i_selected[] = preset.sf)
    apply_preset!(tb_mesh, preset, :mesh)
    apply_preset!(tb_timesteps, preset, :timesteps)
    apply_preset!(tb_quadorder, preset, :quadorder)
    apply_preset!(tb_damp, preset, :damp)
    apply_preset!(tb_sleep, preset, :sleep)

    update_all_funcs()
end

function apply_preset!(textbox::Makie.Textbox, preset, key)
    if key in keys(preset)
        str = getfield(preset, key)
        textbox.displayed_string[] = str
        textbox.stored_string.val = str
    end
end

function apply_bc_preset!(menu::Makie.Menu, preset, key)
    if key in keys(preset)
        i = findfirst(item -> item[2] == getfield(preset, key), menu.options.val)
        menu.i_selected[] = i
    end
end

########################
# %% FEM Model
########################

function update_all_funcs()
    Makie.notify(tb_delta.stored_string)
    Makie.notify(tb_nu.stored_string)
    Makie.notify(tb_gamma.stored_string)
    Makie.notify(tb_beta.stored_string)
    Makie.notify(tb_alpha.stored_string)
    Makie.notify(tb_g.stored_string)
    Makie.notify(tb_bc_left_a.stored_string)
    Makie.notify(tb_bc_left_b.stored_string)
    Makie.notify(tb_bc_left_c.stored_string)
    Makie.notify(tb_bc_right_a.stored_string)
    Makie.notify(tb_bc_right_b.stored_string)
    Makie.notify(tb_bc_right_c.stored_string)
    Makie.notify(tb_dudt_init.stored_string)
    Makie.notify(tb_u_init.stored_string)
end

function getmodel()
    mesh = Mesh(evalparse(tb_mesh), Edge2)

    timesteps = collect(evalparse(tb_timesteps))
    quadorder = parse(Int, tb_quadorder.stored_string[])
    shapefunc = menu_sf.selection[]
    sleeptime = parse(Float64, tb_sleep.stored_string[])
    damping = parse(Float64, tb_damp.stored_string[])

    fields = (
        u = ScalarField{Float64}(mesh, shapefunc, cachelength=1, outputname="Solution field u"),
    )

    pdes = (
        CoefficientFormPDE1D(fields.u; delta, nu, gamma, beta, alpha, g, sleeptime),
    )

    bctype_left = menu_bc_left.selection[]
    bc_left = if bctype_left <: Nothing
        nothing
    elseif bctype_left <: DirichletBC
        DirichletBC(fields.u, :boundary_left, (x, t) -> a_left(NaN, NaN, x, t))
    elseif bctype_left <: GeneralBC
        GeneralBC(fields.u, [:boundary_left],
            (u, ux, x, t) -> a_left(u, ux, x, t),
            (u, ux, x, t) -> b_left(u, ux, x, t),
            (u, ux, x, t) -> c_left(u, ux, x, t),
        )
    end

    bctype_right = menu_bc_right.selection[]
    bc_right = if bctype_right <: Nothing
        nothing
    elseif bctype_right <: DirichletBC
        DirichletBC(fields.u, :boundary_right, (x, t) -> a_right(NaN, NaN, x, t))
    elseif bctype_right <: GeneralBC
        GeneralBC(fields.u, [:boundary_right],
            (u, ux, x, t) -> a_right(u, ux, x, t),
            (u, ux, x, t) -> b_right(u, ux, x, t),
            (u, ux, x, t) -> c_right(u, ux, x, t),
        )
    end

    bcs = Tuple(bc for bc in (bc_left, bc_right) if bc isa AbstractBC)

    fieldhandler = FieldHandler!(fields)
    bchandler = BCHandler(bcs)
    soe = preassemble(pdes, bchandler)

    fields.u.t[] = first(timesteps)
    fields.u.t_cache .= first(timesteps)
    Δt1 = timesteps[2] - timesteps[1]
    interpolate!(fields.u, SpaceTimeFunction((x, t) -> dudt_init(x[1], t)))
    @. fields.u.u_cache[1] = -fields.u.u * Δt1

    interpolate!(fields.u, SpaceTimeFunction((x, t) -> u_init(x[1], t)))
    @. fields.u.u_cache[1] += fields.u.u

    linsolver = LinearSolver(soe, UMFPACKFactorization(), verbose=false)
    nlsolver = NewtonSolver(soe, verbose=true, linsolver=linsolver, damping=damping)
    solver = TransientSolver(nlsolver, timesteps)
    quadrature = GaussQuadrature(quadorder)

    return fields, solver, soe, pdes, fieldhandler, bchandler, quadrature
end

########################
# Execution
########################

menu_preset.i_selected[] = 3
Makie.notify(menu_preset.selection)

# Makie.resize_to_layout!(fig)

display(fig)
