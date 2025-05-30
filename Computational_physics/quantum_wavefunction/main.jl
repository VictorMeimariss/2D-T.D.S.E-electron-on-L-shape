using Plots
using SparseArrays
plotlyjs() # Enable PlotlyJS backend for interactivity
include("functions.jl")
import .Functions

# Victor Emmanuel Meimaris 23/3/25: Created main.jl for project in advanced scientific calc, where I will be solving the time evolving shrodinger equation
# in a 2D L shaped domain, creating my own solves and combining them (bicgstab, domain decomposition and multigrid).

# Input variables
domain = (-1, 1) # In nanometers
time = (0, 1)
domain_min = domain[1]
max_length = 0.5
overlaps = 1 # For domain decomposition
iterations = 1 #3000 #crank nicolson frames

# Wavefunction parameters 
sigma = 0.15
y0 =  - 0.5
x0 = - 0.5
ky =  0.0 # + goes to the negative direction 
kx = -20.0 # - 5.0  

# Create mesh and extract parameters
mesh = Functions.grid(domain, max_length)
coords = mesh[1]
l2g = mesh[2]
noe = mesh[3]
nop = mesh[4]
boundary_nodes = mesh[5]
step_size = mesh[6]
lengthr = mesh[7]

dt = 1e-6 # Time step dt<(dx)^2 for optimal results for ex step_size^2 / 4 


# Domain decomposition matrices 
a = mesh[8]
b = mesh[9]
c = mesh[10]
nx_half = mesh[11]
ny_half = mesh[12]

# Potential function
V_flag = 2 # Potential V flag, 1 = box, 2 = box with circle barrier, 3 = box with circle well
V0 = 15 # Only needed for flags>1 is in eV
x_0 = 0.0 # Circle parameters
y_0 = - 0.25
r_0 = 0.5
V_potential_func = Functions.V_function(V_flag, V0, x_0, y_0, r_0)

# Get matrices
matrices = Functions.fem_matrices(V_potential_func, dt, coords, l2g, noe, boundary_nodes, step_size, lengthr)

println("Number of elements: ", noe)
println("Time step: $dt picoseconds")

# Create wavefunction
psi_0(x, y) = Functions.wavefunction(x, y; x0, y0, sigma, kx, ky)

# Currently commenting out lines i dont need to test my script!

# Save as mp4
# n steps of time
time_domain = time[2] - time[1]
n_steps = Int128(time_domain ÷ dt) + 1

no_frames = 10#17000
frame_inter = 1#100

anim = Functions.animated_solution(coords, nop, psi_0, time, matrices..., no_frames, frame_inter)# or n_steps
mp4(anim, "Animations/Electron/test.mp4", fps=15)
println("Done")


#= Create two plots for testing
psi = Functions.solution(coords, nop, psi_0, time, matrices..., a, b, c, nx_half, ny_half, overlaps, iterations)
return 0

psi_final = abs2.(psi[1])
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
display(p1)=#