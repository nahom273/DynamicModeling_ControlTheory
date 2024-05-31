using PyPlot

# Define control parameters
Λ = 900 # Proportional gain (P)
K = 1/900 #Gain for the corrective action, similar to Derivative gain (D)
w = 1   # Boundary width for the tanh function, smoothing control actions

# Define time parameters
δt = 1e-3
LONG = Int(2e4)
l = LONG - 1

# System dynamics parameters
σₑ = 10
ϱₑ = 100
βₑ = 0.3

# Inverse dynamics parameters for control calculation
σₐ = 9
βₐ = 98
ϱₐ = 0.4

# Nominal trajectory parameters
A1 = 2
ω1 = 0.5
A2 = 3.8
ω2 = 0.7
A3 = 1.5
ω3 = 0.9

# Initialization of state and control arrays
x = zeros(LONG)
y = zeros(LONG)
z = zeros(LONG)

ẋ = zeros(LONG)
ẏ = zeros(LONG)
ż = zeros(LONG)

ẋᴰᵉˢ = zeros(LONG)
ẏᴰᵉˢ = zeros(LONG)
żᴰᵉˢ = zeros(LONG)

xᴺ = zeros(LONG)
yᴺ = zeros(LONG)
zᴺ = zeros(LONG)

ẋᴺ = zeros(LONG)
ẏᴺ = zeros(LONG)
żᴺ = zeros(LONG)

u_x = zeros(LONG)
u_y = zeros(LONG)
u_z = zeros(LONG)

S_x = zeros(LONG)
S_x_p = zeros(LONG)
S_y = zeros(LONG)
S_y_p = zeros(LONG)
S_z = zeros(LONG)
S_z_p = zeros(LONG)

t = zeros(LONG)

e_int_x = 0
e_int_y = 0
e_int_z = 0
e_prev_x = 0
e_prev_y = 0
e_prev_z = 0



# intial conditions
x[1]=0
y[1]=0
z[1]=0

ẋ[1]=0
ẏ[1]=0
ż[1]=0

# Simulation loop
for i = 1:l
    t[i] = i * δt  # Time update
    # Calculate nominal trajectories
    xᴺ[i] = A1 * sin(ω1 * t[i])
    yᴺ[i] = A2 * sin(ω2 * t[i])
    zᴺ[i] = A3 * sin(ω3 * t[i])
    ẋᴺ[i] = A1 * ω1 * cos(ω1 * t[i])
    ẏᴺ[i] = A2 * ω2 * cos(ω2 * t[i])
    żᴺ[i] = A3 * ω3 * cos(ω3 * t[i])
    
    # Calculate current errors
    e_x = xᴺ[i] - x[i]
    e_y = yᴺ[i] - y[i]
    e_z = zᴺ[i] - z[i]
    
    # Integrate error for the Integral component
    e_int_x += δt * e_x
    e_int_y += δt * e_y
    e_int_z += δt * e_z
    
    # Derivative of error for the Derivative component
    e_dot_x = (e_x - e_prev_x) / δt
    e_dot_y = (e_y - e_prev_y) / δt
    e_dot_z = (e_z - e_prev_z) / δt
    e_prev_x = e_x
    e_prev_y = e_y
    e_prev_z = e_z

    # Sliding mode control computation with PID error metrics
    S_x[i] =  e_x + Λ * e_int_x + K * e_dot_x  # Adding PID components
    S_y[i] =  e_y + Λ * e_int_y + K * e_dot_y
    S_z[i] =  e_z + Λ * e_int_z + K * e_dot_z

    
    
    # Compute control inputs using sliding mode dynamics
    ẋᴰᵉˢ[i] = Λ * e_x + ẋᴺ[i] + K * tanh(S_x[i] / w)
    ẏᴰᵉˢ[i] = Λ * e_y + ẏᴺ[i] + K * tanh(S_y[i] / w)
    żᴰᵉˢ[i] = Λ * e_z + żᴺ[i] + K * tanh(S_z[i] / w)

   # Dynamics inverse model
    u_x[i] = ẋᴰᵉˢ[i] - σₐ * y[i]
    u_y[i] = ẏᴰᵉˢ[i] + βₐ * x[i] - y[i] * z[i]
    u_z[i] = żᴰᵉˢ[i] - ϱₐ + y[i]^2

 # Dynamics Update system states using the system model
    ẋ[i] = σₑ * y[i] + u_x[i]
    ẏ[i] = -βₑ * x[i] + y[i] * z[i] + u_y[i]
    ż[i] = ϱₑ - y[i]^2 + u_z[i]
    
    x[i+1] = x[i] + δt * ẋ[i]
    y[i+1] = y[i] + δt * ẏ[i]
    z[i+1] = z[i] + δt * ż[i]
end

# Plot the results: Nominal and Realized Trajectories, Tracking Errors, Phase Spaces, Control Signals
using PyPlot

# Plotting nominal and realized trajectories
figure(1)
grid(true)
title("Nominal and Realized Trajectory")
plot(t[1:l], xᴺ[1:l], label="Nominal x")
plot(t[1:l], yᴺ[1:l], label="Nominal y")
plot(t[1:l], zᴺ[1:l], label="Nominal z")
plot(t[1:l], x[1:l], "--", label="Realized x")
plot(t[1:l], y[1:l], "--", label="Realized y")
plot(t[1:l], z[1:l], "--", label="Realized z")
legend()
savefig("Nominal and Realized Trajectory.pdf")

# Plotting tracking errors for x, y, z
figure(2)
grid(true)
title("Tracking Error")
plot(t[1:l], xᴺ[1:l] - x[1:l], label="Error x")
plot(t[1:l], yᴺ[1:l] - y[1:l], label="Error y")
plot(t[1:l], zᴺ[1:l] - z[1:l], label="Error z")
legend()
savefig("Tracking Error.pdf")

# Plotting phase space trajectories
figure(3,figsize=(15, 5))
grid(true)
title("Phase Space")
plot(ẋᴺ[1:l], xᴺ[1:l], label="Nominal Phase x")
plot(ẏᴺ[1:l], yᴺ[1:l], label="Nominal Phase y")
plot(żᴺ[1:l], zᴺ[1:l], label="Nominal Phase z")
plot(ẋ[1:l], x[1:l], "--", label="Realized Phase x")
plot(ẏ[1:l], y[1:l], "--", label="Realized Phase y")
plot(ż[1:l], z[1:l], "--", label="Realized Phase z")
legend()
savefig("Phase Space.pdf")

# Plotting control signals
figure(4)
grid(true)
title("Control Signals")
plot(t[1:l], u_x[1:l], label="Control u_x")
plot(t[1:l], u_y[1:l], label="Control u_y")
plot(t[1:l], u_z[1:l], label="Control u_z")
legend()
savefig("Control Signals.pdf" ) # add this into update it later to latex

show()
