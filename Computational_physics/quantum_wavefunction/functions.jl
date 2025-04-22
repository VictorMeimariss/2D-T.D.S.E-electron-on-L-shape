#= Victor Emmanuel Meimaris 23/3/25: Creating Module with all functions used on main.jl, I will be
explaining the numbers in my report. =#
module Functions
using SparseArrays, LinearAlgebra, Arpack
using Plots
using SparseArrays
plotlyjs() # Enable PlotlyJS backend for interactivity

# Constants
const h_bar = 1#1.054571817 * 10^(-34)
const m = 1#9.1093837 * 10^(-31)
const gamma = (h_bar^2)/(2*m)

# Creating Coordinate/Mesh function, returns grid ,l2g matrix and number of elements / nodes.
function grid(domain::Tuple{Real, Real}, max_length::Float64) # num_of_squares_x::UInt16, num_of_squares_y::UInt16

    # Define domain length in any axis since x and y have the same domain
    domain_length = domain[2] - domain[1]

    # Filtering input so that it's valid 

    # This division has to be a multiple of 2 in order for FEM to work so we modify max_length accordingly
    temp = domain_length / max_length
    if temp < 2
        max_length = domain_length / 2
    elseif temp % 2 != 0
        temp = round(Int, temp)  # Round to the nearest integer
        if temp % 2 != 0
            temp += 1
        elseif temp == 2 # Added this for occasions where temp rounds to 2 so it needs +2
            temp += 2
        end
        max_length = domain_length / temp
    end
    temp = round(Int, temp) # Case for if temp = 2.0

    # Length of xg and yg
    lengthr = Int64(domain_length / max_length + 1)

    # Deriving number of squares based on length provided
    num_of_squares_y = temp
    num_of_squares_x = temp # This could be different to y if multiple of 2

    # Number of elements and nodes
    noe = Int((num_of_squares_y * num_of_squares_x) * 0.75) # Number of elements
    nop = num_of_squares_x + num_of_squares_y + noe + 1 # Number of nodes

    # a->top left square
    a = reshape(1:((num_of_squares_x ÷ 2) + 1) * ((num_of_squares_y ÷ 2) + 1), (num_of_squares_x ÷ 2) + 1, (num_of_squares_y ÷ 2) + 1)
    a = transpose(a)
    
    # b->bottom big rectangle,b and a overlap on a's last line, for convinience in when creating the l2g matrix
    b = reshape(a[end,1]:nop, num_of_squares_x + 1, (num_of_squares_y ÷ 2) + 1)
    b = transpose(b)

    # Keep boundary nodes for boundary conditions in the solution later without any order 
    boundary_nodes = Set{Int}() # Initialise empty set for boundary nodes for fast look ups to detect number of boundary nodes in each element
    union!(boundary_nodes, @view a[1, :])
    union!(boundary_nodes, @view a[:, 1])
    union!(boundary_nodes, @view a[:, end])
    union!(boundary_nodes, @view b[:, 1])
    union!(boundary_nodes, @view b[end, :])
    union!(boundary_nodes, @view b[:, end])
    union!(boundary_nodes, @view b[1,((num_of_squares_x ÷ 2) + 1) :end])

    # Generate coordinate range
    x = range(domain[1], domain[2], num_of_squares_x + 1)# need boundary nodes for dirichlet, elements with two sides with boundary nodes and elem with 1 side
    y = range(domain[2], domain[1], num_of_squares_y + 1)#if boundary nodes in l2g ==2 then 1 side elseif ==3 2 sides else no sides create flag or smthing to not do ifs all the time

    # Generate coordinates
    xc = [xi for xi in x, _ in y]  # Each row is xi, number of rows = length of y
    yc = [yi for _ in x, yi in y]  # Each column is yi, number of collumns = length of x

    # Setting coordinates on a grid
    coords = zeros(nop, 2)
    temp = 1
    index = 1
    for i = 1:(num_of_squares_y ÷ 2)
        @inbounds for j = 1:(num_of_squares_x ÷ 2 + 1)
            coords[index, 1] = xc[temp]
            coords[index, 2] = yc[temp]
            temp += 1
            index += 1
        end
        temp += (num_of_squares_y ÷ 2) # I am not adding +1 to temp or index because in last iteration of j loop it was already done
    end
    
    coords[index:end, 1] = xc[temp:length(xc)]
    coords[index:end, 2] = yc[temp:length(yc)]

    # Define local to global map " l2g "
    l2g = zeros(Int, noe, 4) # Each element has 4 nodes
    index = 1
    # For top square
    for i = 1:num_of_squares_x ÷ 2
        @inbounds for j = 1:num_of_squares_y ÷ 2
            l2g[index, :] = [a[i + 1, j], a[i + 1, j + 1], a[i, j + 1], a[i, j]]
            index +=1
        end
    end
    # For bottom rectangle
    for i = 1:(num_of_squares_y ÷ 2)
        @inbounds for j = 1: num_of_squares_x
            l2g[index, :] = [b[i + 1, j], b[i + 1, j + 1], b[i, j + 1], b[i, j]]
            index += 1
        end
    end
    return coords, l2g ,noe, nop, boundary_nodes, max_length, lengthr
end

#= Creating function which will be returning the matrices A and B containing the Hamiltonian and Mass matrices using F.E.M
    A and B are the LHS and RHS matrices used in the crank nicolson solution to provide wavefunction solution
=#
function fem_matrices(V_potential_func, dt, coords::Matrix{Float64}, l2g::Matrix{Int64}, noe::Int64, nop::Int64, boundary_nodes::Set, step_size, lengthr)# step _size is here because its needed for plotting
    
    # For this script I will assume Neuman and robin coefficients to be zero and then might integrate them for certain problems like scattering
    # I will also assume Dirichlet boundary conditions where the probability ψ is zero in the boundaries(Potential V = infinite in the boundaries), as in this general function I want to model
    # a closed system. If an open system would be simulated, the matrix p would be needed as well, whichs contains the neuman and robin coefficients

    # Shape/basis functions for integration
    N = [
        (ksi,eta) -> 1/4 * (1 - ksi) * (1 - eta),# N1
        (ksi,eta) -> 1/4 * (1 + ksi) * (1 - eta),# N2
        (ksi,eta) -> 1/4 * (1 + ksi) * (1 + eta),# N3
        (ksi,eta) -> 1/4 * (1 - ksi) * (1 + eta)# N4
    ]

    # Assemble Hamiltonian and Mass matrix indices (triplets for sparse matrix) H = K + V matrices
    ia = zeros(noe * 16) # Row index, noe * 16 because the two for loops for each square = 16 iterations for each element --> noe * 16
    ja = zeros(noe * 16) # Column index
    va_H = zeros(noe * 16) # Hamiltonian matrix value index
    va_M = zeros(noe * 16) # Mass matrix value index

    # Iterate over elements to find global stiffness and force matrices
    index = 1
    for e = 1:noe

        xe = coords[l2g[e, :], 1]
        ye = coords[l2g[e, :], 2]

        # Local to globall coordinates
        x(ksi, eta) = N[1](ksi, eta) * xe[1] + N[2](ksi, eta) * xe[2] + N[3](ksi, eta) * xe[3] + N[4](ksi, eta) * xe[4]
        y(ksi, eta) = N[1](ksi, eta) * ye[1] + N[2](ksi, eta) * ye[2] + N[3](ksi, eta) * ye[3] + N[4](ksi, eta) * ye[4]

        # Creating jacobian matrix for element e
        J11(eta) = 0.25 * (- (1 - eta) * xe[1] + (1 - eta) * xe[2] + (1 + eta) * xe[3] - (1 + eta) * xe[4])
        J21(ksi) = 0.25 * (- (1 - ksi) * xe[1] - (1 + ksi) * xe[2] + (1 + ksi) * xe[3] + (1 - ksi) * xe[4])
        J12(eta) = 0.25 * (- (1 - eta) * ye[1] + (1 - eta) * ye[2] + (1 + eta) * ye[3] - (1 + eta) * ye[4])
        J22(ksi) = 0.25 * (- (1 - ksi) * ye[1] - (1 + ksi) * ye[2] + (1 + ksi) * ye[3] + (1 - ksi) * ye[4])
        detJ(ksi, eta) = J11(eta) * J22(ksi) - J12(eta) * J21(ksi)
        

         # Define θNx and θNy as arrays of functions
         θNx = [
            (ksi, eta) -> 0.25 * (-J22(ksi) * (1 - eta) + J12(eta) * (1 - ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (J22(ksi) * (1 - eta) + J12(eta) * (1 + ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (J22(ksi) * (eta + 1) - J12(eta) * (ksi + 1)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (-J22(ksi) * (eta + 1) - J12(eta) * (1 - ksi)) / detJ(ksi, eta)
        ]
        θNy = [
            (ksi, eta) -> 0.25 * (J21(ksi) * (1 - eta) - J11(eta) * (1 - ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (-J21(ksi) * (1 - eta) - J11(eta) * (1 + ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (-J21(ksi) * (1 + eta) + J11(eta) * (1 + ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (J21(ksi) * (1 + eta) + J11(eta) * (1 - ksi)) / detJ(ksi, eta)
        ]
        # Local Hamiltonian and Mass matrices
        KV = zeros(4, 4)
        L = zeros(4, 4) # Local Mass matrix
        
        # Precompute integrands for less allocations
        # Local potential matrix V integrals
        integrands_V = [
            (ksi, eta) -> N[i](ksi, eta) * N[j](ksi, eta) * V_potential_func(x(ksi, eta), y(ksi, eta)) * abs(detJ(ksi, eta))
            for i in 1:4, j in 1:4
        ]

        # Compute local stiffness matrix K integrals with αx and αy equal to 1 for the Shrodinger equation
        integrands_K = [
            (ksi, eta) -> gamma * (θNx[i](ksi, eta) * θNx[j](ksi, eta) + θNy[i](ksi, eta) * θNy[j](ksi, eta)) * abs(detJ(ksi, eta))
            for i in 1:4, j in 1:4
        ]

        # Compute local Mass matrix integrals
        integrands_L = [
            (ksi, eta) -> N[i](ksi, eta) * N[j](ksi, eta) * abs(detJ(ksi, eta))
            for i in 1:4, j in 1:4
        ]
        # Add results to local Hamiltonian and Mass matrix 
        for i = 1:4
            @inbounds for j = 1:4
                KV[i, j] = gauss_quad_2D(integrands_K[i, j], 3) + gauss_quad_2D(integrands_V[i, j], 3)
                L[i, j] = gauss_quad_2D(integrands_L[i, j], 3)
            end
        end
        
        # Assemble global stiffness and Mass matrices indices from the local element matrices 
        for i = 1:4
            @inbounds for j = 1:4
                ia[index] = l2g[e, i] # global node index at row i of local matrices
                ja[index] = l2g[e, j] # global node index at collumn j
                va_H[index] = KV[i, j] # global value index at i, j of local Hamiltonian matrix
                va_M[index] = L[i, j] # global value index at i, j of local Mass matrix
                index = index + 1;
            end
        end
    end
    # Assembling global Hamiltonian matrix = global stiff+ global potential, and Mass matrices using indices from loop
    H = sparse(ia, ja, va_H)
    M = sparse(ia, ja, va_M)
    
    # Enforce dirichlet boundary conditions
    boundary_nodes = collect(boundary_nodes)  # Convert Set to Vector
    for node in boundary_nodes
        H[node, :] .= 0  # Zero row
        H[:, node] .= 0  # Zero column
        M[node, :] .= 0
        M[:, node] .= 0
        H[node, node] = 1  # Set diagonal
        M[node, node] = 1
    end
    #= Tests
    println("M symmetric?",norm(M-M')<1e-10)
    println("H symmetric?",norm(H-H')<1e-10)
    
    # Test H for positive semi-definiteness
    try
        lambda_min = eigs(H, nev=1, which=:SR)[1][1]
        lambda_max = eigs(H, nev=1, which=:LR)[1][1]
        if real(lambda_min) >= -1e-10
            println("H is positive semi-definite, smallest eigenvalue: ", lambda_min)
            println("H highest eigenvalue: ",lambda_max)
        else
            println("H is not positive semi-definite, smallest eigenvalue: ", lambda_min)
        end
    catch e
        println("Error computing eigenvalues of H: ", e)
    end
    try
        chol_M = cholesky(M)
        println("M is positive definite")
    catch e
        println("M is not positive definite: ", e)
    end=#
    # Multiply H by i * dt / (2*h_bar) to include everything in H and hold previous value in tempo for latter use
    tempo = H
    H = H * im * dt / (2 * h_bar)

    # Matrices A and B for RHS and LHS
    A = H + M
    B = M - H
    
    # Now our system looks like this: Aψ(n+1) = Βψn, almost ready for the solution with crank nicolson
    return A, B, lengthr, dt, M, boundary_nodes, tempo
end

# Creating wavefunction equation to solve with function solution
function wavefunction(x, y ; x0, y0, sigma_x, sigma_y, kx_unscaled::Float64, ky_unscaled::Float64, step_size)
    kx = kx_unscaled * pi / 2#* ((2 * pi) / (1 / step_size))
    ky = ky_unscaled * pi / 2#* ((2 * pi) / (1 / step_size))
    gaussian = (1 / (pi * sigma_x * sigma_y)) * exp(-((x - x0)^2 / (2 * sigma_x^2) + (y - y0)^2 / (2 * sigma_y^2)))
    plane_wave = exp(im * (kx * x + ky * y))
    return gaussian * plane_wave
end
# Creating function for getting solution from Ax=b where A is A,x is ψ(n+1) and b = B*ψn
function solution(coords, nop, psi_zero, time, A, B, lengthr, dt, M, boundary_nodes, tempo)
    
    # n steps of time
    time_domain = time[2] - time[1]
    n_steps = Int128(time_domain ÷ dt )+ 1

    # Initialize psi_0 as a 1D vector
    psi_0 = Vector{ComplexF64}(undef, nop)

    # Extracting data from psi_zero function to psi_0 vector on our coordinate system
    for k in 1:nop
        psi_0[k] = psi_zero(coords[k,1], coords[k,2])
    end

    # Enforce Dirichlet boundary conditions (ψ = 0 at boundaries)
    psi_0[collect(boundary_nodes)] .= 0

    norm = sqrt(real(psi_0' * M * psi_0))
    #println("Unnormalized max |ψ₀|: ", maximum(abs.(psi_0)))
    #println("Unnormalized norm: ", norm)
    
    psi_0 = psi_0 / norm
    temp = psi_0
    #println("Normalized max |ψ₀|: ", maximum(abs.(psi_0)))
    #println("Normalized norm check: ", sqrt(real(psi_0' * M * psi_0)))

    # Calculate initial energy
    #E_initial = real(psi_0' * tempo * psi_0)
    #println("Initial energy: ", E_initial)

    # Time stepping setup
    A_LU = lu(A)
    psi = similar(psi_0)
    # Time stepping loop
    for n = 1:10000#n_steps

        # Solve psi
        psi = A_LU \ (B * psi_0)
        
        # Asigning value to psi_0
        psi_0 = psi

    end
    
    #= Final calculations
    final_norm = sqrt(real(psi_0' * M * psi_0))
    final_E = real(psi_0' * tempo * psi_0) / final_norm
    
    println("\nFinal results:")
    println("Norm: $final_norm")
    println("Max |ψ|: ", maximum(abs.(psi_0)))
    println("Energy: $final_E")=#

    return psi, temp # return psi and the initial state
end
function animated_solution(step_size, coords, nop, psi_zero, time, A, B, lengthr, dt, M, boundary_nodes, tempo)
    # Initialize
    psi_0 = Vector{ComplexF64}(undef, nop)

    # Extracting data from psi_zero function to psi_0 vector on our coordinate system
    for k in 1:nop
        psi_0[k] = psi_zero(coords[k,1], coords[k,2])
    end

    # Enforce Dirichlet boundary conditions (ψ = 0 at boundaries)
    psi_0[collect(boundary_nodes)] .= 0

    norm = sqrt(real(psi_0' * M * psi_0))
    #println("Unnormalized max |ψ₀|: ", maximum(abs.(psi_0)))
    #println("Unnormalized norm: ", norm)
    
    psi_0 = psi_0 / norm
    
    # Precompute
    A_LU = lu(A)
    frames = []
    xg = yg = range(-1, 1, step=step_size)

    # Calculate global maximum for consistent scaling
    global_max = maximum(abs2.(psi_0)) * 1.5 + 2
    
    # Time evolution
    for n in 1:120000 # Reduce frames for quicker animation
        psi = A_LU \ (B * psi_0)
        
        # Capture every 200th frame
        if n % 200 == 0
            Z = fill(NaN, lengthr, lengthr)
            for k in 1:nop
                i, j = coord_to_index(coords[k,:], step_size)
                Z[i,j] = abs2(psi[k])
            end
            push!(frames, surface(
                xg, yg, Z,
                title = "Wavefunction",
                colorbar = true,
                colorscale = "Viridis",
                showscale = true,
                
                # Critical stabilization parameters:
                zlim = (0, global_max),          # Fix z-axis range
                xlim = (minimum(xg), maximum(xg)), # Fix x-axis
                ylim = (minimum(yg), maximum(yg)), # Fix y-axis
                clims = (0, global_max),         # Lock color scale
                # Backend-agnostic stabilization
                aspect_ratio = :equal,      # Lock proportions
                legend = false,             # Remove legend
                grid = false,              # Remove grid lines (optional)
                
                camera = (30, 25, 1.0),       # Fixed view angle
                show_boundingbox = true, # Force consistent box rendering
                
                # These prevent any automatic adjustments:
                autoscale = :none,
                normalize = false,
                ))
        end
        psi_0 = psi
    end
    return frames
end
# Helper function for animation
function coord_to_index(coord, step_size)
    i = round(Int, (coord[1] + 1) / step_size) + 1
    j = round(Int, (coord[2] + 1) / step_size) + 1
    return (i, j)
end
# Creating 2 dimensional gauss_quadrature function to use to solve the integrals needed in FEM integrating from -1 to 1
function gauss_quad_2D(funct, n)# omada, mission planner, api, 
    if n == 2
        ksi = (-0.5773502692, 0.5773502692)
        eta = (-0.5773502692, 0.5773502692)
        w = ones(2)
    elseif n == 3
        ksi = (0, -0.7745966692, 0.7745966692)
        eta = (0, -0.7745966692, 0.7745966692)
        w = (0.8888888889, 0.5555555556, 0.5555555556)
    else
        error("Only n = 2 or n = 3 are supported")
    end
    result = 0.0
    for i = 1:n
        @inbounds for j = 1:n
            result += w[i] * w[j] * funct(ksi[i], eta[j])
        end
    end
    return result
end
# Creating function for potential V and its flags
function V_function(V_flag::Int64, V0)
    if V_flag == 1 # Infinite potential well with bottom at V0
        return (x,y) -> V0
    elseif V_flag == 2 # Barrier of potential V0 at area around x0 and y0
        return (x, y) -> V0
    elseif V_flag == 3 # Well of potential V0
        return (x, y) -> -V0
    end
end
end # Module end