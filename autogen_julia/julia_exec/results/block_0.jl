using Plots
using DifferentialEquations

# Pendulum parameters
pendulum_length = 1.0      # Length of the pendulum stick
gravity = 9.81             # Acceleration due to gravity
damping_coefficient = 0.1  # Damping coefficient

# Initial conditions
initial_angle = π / 4      # Initial angle in radians
initial_angular_velocity = 0.0 # Initial angular velocity

# Define the pendulum equations with damping
function pendulum_equations!(du, u, params, t)
    theta, omega = u             # θ is angle, ω is angular velocity
    du[1] = omega                # dθ/dt = ω
    du[2] = -(params[1] / params[2]) * sin(theta) - params[3] * omega  # dω/dt
end

# Setup and solve the differential equation
g = gravity
L = pendulum_length
b = damping_coefficient
params = [g, L, b]
u0 = [initial_angle, initial_angular_velocity]  # Initial state
tspan = (0.0, 10.0)  # Time span for simulation

prob = ODEProblem(pendulum_equations!, u0, tspan, params)
sol = solve(prob, Tsit5())

# Create and save a gif of the pendulum motion
anim = @animate for i in 1:length(sol)
    angle = sol.u[i][1]  # Use the correct solution to get current angle
    x = pendulum_length * sin(angle)
    y = -pendulum_length * cos(angle)

    plot([0, x], [0, y], line = (2, :blue), size=(400, 400),
         xlims=(-L, L), ylims=(-L, L), aspect_ratio=:equal)
    scatter!([x], [y], color = :red, markersize = 6)
end

gif(anim, "pendulum_simulation.gif", fps=20)  # Directly pass the filename and fps