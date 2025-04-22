using Plots
using SparseArrays
#gr()
plotlyjs() # Enable PlotlyJS backend for interactivity
#plotly()
include("functions.jl")
import .Functions

# Input variables
domain = (-1, 1)
time = (0, 1)
dt = 1e-5 # Time step dt<(dx)^2 for results
domain_min = domain[1]
max_length = 0.05
V_flag = 1 # Potential V flag

# Wavefunction parameters
sigma_x = 0.4
sigma_y = 0.4
x0 = 0.5
y0 = -0.5
kx_unscaled = -3.0
ky_unscaled = 0.0

# Create mesh and extract parameters
mesh = Functions.grid(domain, max_length)
coords = mesh[1]
l2g = mesh[2]
noe = mesh[3]
nop = mesh[4]
boundary_nodes = mesh[5]
step_size = mesh[6]
lengthr = mesh[7]

# Potential function
V_potential_func = Functions.V_function(V_flag, 0.0)
# Create grid ranges
xg = range(domain_min, -domain_min, step=step_size)
yg = range(domain_min, -domain_min, step=step_size)

# Get matrices
matrices = Functions.fem_matrices(V_potential_func, dt, mesh...)

println("Number of elements: ", noe)
# Initialize Z and F matrices
Z = fill(NaN, lengthr, lengthr) # Fill with NaN so that no space is occupied
F = fill(NaN, lengthr, lengthr)
# Define wavefunction and its real part
psi_0(x, y) = Functions.wavefunction(x, y; x0, y0, sigma_x, sigma_y, kx_unscaled, ky_unscaled, step_size)


#=psi = Functions.solution(coords, nop, psi_0, time, matrices...)
psi_final = abs2.(psi[1])#real(psi[1])#
psi_initial = abs2.(psi[2])#real(psi[2])#

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
        xlabel="x", ylabel="y", zlabel="|ψ0|^2",
        camera=(30, 40),
        colorbar=true,
        colorscale="Viridis",
        showscale=true)
display(p2)
p1 = surface(xg, yg, Z, 
        title="Final Wavefunction |ψ0|^2", 
        xlabel="x", ylabel="y", zlabel="|ψf|^2",
        camera=(30, 40),
        colorbar=true,
        colorscale="Viridis",
        showscale=true)
display(p1)=#

# Generate animation
println("Creating animation...")
frames = Functions.animated_solution(step_size, coords, nop, psi_0, time, matrices...)
lengthf = length(frames)
# Save as mp4
anim = @animate for (i, frame) in enumerate(frames)
        plot(frame, camera = (30, 25, 1.0))  # Plot the current frame
        println("Progress: $i/$lengthf")  # Show current/total
end
mp4(anim, "wavefunction_evolution_8.mp4", fps=30)
println("Done")