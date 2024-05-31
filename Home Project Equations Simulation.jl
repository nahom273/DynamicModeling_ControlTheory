using PyPlot

#Parameters (approximations)
a_appr = 9
b_appr = 98
c_appr = 0.4

# Parameters (this are exact parameters. Duirng the class we tired  with exact parameters)
a = 10
b = 100
c = 0.3



# Simulation settings
LONG = Int(2e4)  # Simulation length
l = LONG - 1
h = 1e-3  # Time step

# Initial conditions
x = zeros(LONG)
y = zeros(LONG)
z = zeros(LONG)

x[1] = 0.1
y[1] = 0.1
z[1] = 0.1

# Simulation using Euler's method
for i in 1:l
    ẋ = a*y[i]  # Modified differential equation remove input control as per instuction 
    ẏ = -c*x[i] + y[i]*z[i]  # Modified differential equation
    ż = b - y[i]^2  # Modified differential equation

    x[i+1] = x[i] + ẋ*h
    y[i+1] = y[i] + ẏ*h
    z[i+1] = z[i] + ż*h
end

# Plotting the results
figure(figsize=(10, 8))
subplot(2, 2, 1)
plot3D(x[1:l], y[1:l], z[1:l])
title("3D Plot  ")

subplot(2, 2, 2)
plot(x[1:l])
title("Time plot for x")

subplot(2, 2, 3)
plot(y[1:l])
title("Time plot for y")

subplot(2, 2, 4)
plot(z[1:l])
title("Time plot for z")

tight_layout()
show()

