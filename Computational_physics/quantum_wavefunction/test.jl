using Plots
using SparseArrays
plotlyjs() # Enable PlotlyJS backend for interactivity
include("functions.jl")
import .Functions

# Input variables
domain = (-1, 1) # In nanometers
time = (0, 1)
domain_min = domain[1]
max_length = 0.025

# Wavefunction parameters
sigma = 0.15
y0 =  - 0.5
x0 = - 0.5
ky =  0.0 # + goes to the negative direction 
kx = -20.0 #- 5.0

# Create mesh and extract parameters
mesh = Functions.grid(domain, max_length)
coords = mesh[1]
l2g = mesh[2]
noe = mesh[3]
nop = mesh[4]
boundary_nodes = mesh[5]
step_size = mesh[6]
lengthr = mesh[7]
dt = 1e-6#step_size^2 / 4 # Time step dt<(dx)^2 for results

# Potential function
V_flag = 3 # Potential V flag
V0 = 7# Only needed for flags>1 7eV
x_0 = 0.0
y_0 = - 0.5
r_0 = 0.5
V_potential_func = Functions.V_function(V_flag, V0, x_0, y_0, r_0)

# Create grid ranges
xg = range(domain_min, -domain_min, step=step_size)
yg = range(domain_min, -domain_min, step=step_size)

# Get matrices
matrices = Functions.fem_matrices(V_potential_func, dt, mesh...)

println("Number of elements: ", noe)
println("Time step: ", dt) 
# Initialize Z and F matrices
Z = fill(NaN, lengthr, lengthr) # Fill with NaN so that no space is occupied
F = fill(NaN, lengthr, lengthr)
# Define wavefunction and its real part
psi_0(x, y) = Functions.wavefunction(x, y; x0, y0, sigma, kx, ky)

# Currently commenting out lines i dont need to test my script!


# Save as mp4
anim = Functions.animated_solution(coords, nop, psi_0, time, matrices..., 17000, 100)
mp4(anim, "Animations/Electron/electron_x_potential well_7eV.mp4", fps=15)#infinite_potential_well_2
println("Done")


#= Create two plots for testing
psi = Functions.solution(coords, nop, psi_0, time, matrices...)
#=psi_final = abs2.(psi[1])
psi_initial = abs2.(psi[2])

# Define grid ranges and Z
xg = -1:step_size:1
yg = -1:step_size:1
Z = fill(NaN, lengthr, lengthr) # Fill with NaN so that no space is occupied
F = fill(NaN, lengthr, lengthr) 
# Map solution vector to grid
for k in 1:nop # Length of solution vector
    x, y = coords[k, 1], coords[k, 2]
    i = round(Int, (x - domain_min) / step_size) + 1  # +1 because Julia is 1-based
    j = round(Int, (y - domain_min) / step_size) + 1
    Z[i, j] = psi_final[k]
    F[i, j] = psi_initial[k]
end

# Create plots
p2 = surface(xg, yg, F, 
        title="Initial Wavefunction |ψ0|^2", 
        xlabel="y", ylabel="x", zlabel="|ψ0|^2",
        camera=(30, 40),
        colorbar=true,
        colorscale="Viridis",
        showscale=true)
display(p2)
p1 = surface(xg, yg, Z, 
        title="Final Wavefunction |ψ0|^2", 
        xlabel="y", ylabel="x", zlabel="|ψf|^2",
        camera=(30, 40),
        colorbar=true,
        colorscale="Viridis",
        showscale=true)
display(p1)=#=#