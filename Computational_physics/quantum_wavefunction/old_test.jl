using Plots
using SparseArrays
plotlyjs() # Enable PlotlyJS backend for interactivity
include("old_functions_test.jl")
import .Old_Functions

# Input variables
domain = (-1, 1)
domain_min = domain[1] # Used later for grid creation
max_length = 0.02
V_flag = 1 # Potential V flag for different potentials of Shrodinger equation

mesh = Old_Functions.grid(domain, max_length)

# Extracting parameters for plotting
coords = mesh[1]
l2g = mesh[2]
noe = mesh[3]
nop = mesh[4]
step_size = mesh[6]

# Right hand side equation test
g(x, y) = -2 * pi^2 .* sin(pi * x) .* sin(pi * y)

# Creating matrix equation
matrix_equation = Old_Functions.matrix_equation(g, mesh..., V_flag)
# Solution vector
u = Old_Functions.solution(matrix_equation...)

# Define grid ranges and Z
xg = -1:step_size:1
yg = -1:step_size:1
Z = fill(NaN, length(xg), length(yg)) # Fill with NaN so that no space is occupied

# Map solution vector to grid
for k in 1:nop # Length of solution vector
    x, y = coords[k, 1], coords[k, 2]
    i = round(Int, (x - domain_min) / step_size) + 1  # +1 because Julia is 1-based
    j = round(Int, (y - domain_min) / step_size) + 1
    Z[i, j] = u[k]
end

println("Number of elements: ", noe)
println("Error at max: ", 1 - maximum(u))

# Plot with interactive surface
p1 = surface(xg, yg, Z, 
        title="Solution u on L-shaped Domain", 
        xlabel="x", ylabel="y", zlabel="u",
        camera=(30, 40), # Adjust camera angle for better view
        colorbar=true,
        colorscale="Viridis", # Colormap
        showscale=true)
display(p1)# Display allows for more than one plots, I will be plotting 3
# First we will solve the time independant shrodinger equation and then integrate time later with finite differences