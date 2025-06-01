#= Victor Emmanuel Meimaris 23/3/25: Creating Module with all functions used on main.jl, I will be
explaining the numbers in my report. =#
module Functions
using SparseArrays, LinearAlgebra, Arpack, Plots, SparseArrays, IncompleteLU

plotlyjs() # Enable PlotlyJS backend for interactivity
#gr() # Backend for animations change to this for stable animations as well as low memory images when testing

# Electron constants for my domain in nanometers and femtoseconds
const h_bar = 0.65821220e-3 # h_bar actual value in eV * ps  
const m = 9.1093837e-31 * (10e25 / 1.602117662 )# mass of electron in kg is 9.1093837e-31 but we want in (eV * ps^2) / nm^2
const gamma = (h_bar^2)/(2*m) # in eV * nm^2

# Creating Coordinate/Mesh function, returns grid ,l2g Matrix and number of elements / nodes.
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
    nop = num_of_squares_x + num_of_squares_y + noe + 1# Number of nodes

    nx_half = (num_of_squares_x ÷ 2) + 1
    ny_half = (num_of_squares_y ÷ 2) + 1
    # a->top left square
    a = reshape(1:nx_half * ny_half, nx_half, ny_half)
    a = transpose(a)
    
    # b->bottom big rectnagle
    b = reshape(a[end,1]:nop, num_of_squares_x + 1, (num_of_squares_y ÷ 2) + 1)
    b = transpose(b)

    # c->bottom right square
    c = b[:, nx_half:end]
    # b -> becomes bottom left square
    b = b[:, 1:nx_half]

    # Keep boundary nodes for boundary conditions in the solution later without any order 
    boundary_nodes = Set{Int}() # Initialise empty set for boundary nodes for fast look ups to detect number of boundary nodes in each element
    union!(boundary_nodes, @view a[1, :])
    union!(boundary_nodes, @view a[:, 1])
    union!(boundary_nodes, @view a[:, end])
    union!(boundary_nodes, @view b[:, 1])
    union!(boundary_nodes, @view b[end, :])
    union!(boundary_nodes, @view c[1, :])
    union!(boundary_nodes, @view c[end, :])
    union!(boundary_nodes, @view c[:, end])

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
    # For bottom left square
    for i = 1:(num_of_squares_y ÷ 2)
        @inbounds for j = 1:(num_of_squares_x  ÷ 2)
            l2g[index, :] = [b[i + 1, j], b[i + 1, j + 1], b[i, j + 1], b[i, j]]
            index += 1
        end
    end
    # For bottom right square overlapping elements
    for i = 1:(num_of_squares_y ÷ 2)
        @inbounds for j = 1:(num_of_squares_x  ÷ 2)
            l2g[index, :] = [c[i + 1, j], c[i + 1, j + 1], c[i, j + 1], c[i, j]]
            index += 1
        end
    end
    return coords, l2g, noe, nop, boundary_nodes, max_length, lengthr, a, b, c, nx_half, ny_half
end

#= Creating function which will be returning the matrices A and B containing the Hamiltonian and Mass matrices using F.E.M
    A and B are the LHS and RHS matrices used in the crank nicolson solution to provide wavefunction solution
=#
function fem_matrices(V_potential_func, dt, coords::Matrix{Float64}, l2g::Matrix{Int64}, noe::Int64, boundary_nodes::Set, step_size, lengthr)# step _size is here because its needed for plotting
    
    # For this script I will assume Neuman and robin coefficients to be zero and then might integrate them for certain problems like scattering
    # I will also assume Dirichlet boundary conditions where the probability ψ is zero in the boundaries(Potential V = infinite in the boundaries), as in this general function I want to model
    # a closed system. If an open system would be simulated, the Matrix p would be needed as well, whichs contains the neuman and robin coefficients
    println("Creating F.E.M matrices...")

    time0 = Base.time()
    # Shape/basis functions for integration
    N = [
        (ksi,eta) -> 1/4 * (1 - ksi) * (1 - eta),# N1
        (ksi,eta) -> 1/4 * (1 + ksi) * (1 - eta),# N2
        (ksi,eta) -> 1/4 * (1 + ksi) * (1 + eta),# N3
        (ksi,eta) -> 1/4 * (1 - ksi) * (1 + eta)# N4
    ]

    # Assemble Hamiltonian and Mass Matrix indices (triplets for sparse Matrix) H = K + V matrices
    ia = zeros(noe * 16) # Row index, noe * 16 because the two for loops for each square = 16 iterations for each element --> noe * 16
    ja = zeros(noe * 16) # Column index
    va_H = zeros(noe * 16) # Hamiltonian Matrix value index
    va_M = zeros(noe * 16) # Mass Matrix value index

    # Iterate over elements to find global stiffness and force matrices
    index = 1
    for e = 1:noe

        xe = coords[l2g[e, :], 1]
        ye = coords[l2g[e, :], 2]

        # Local to globall coordinates
        x(ksi, eta) = N[1](ksi, eta) * xe[1] + N[2](ksi, eta) * xe[2] + N[3](ksi, eta) * xe[3] + N[4](ksi, eta) * xe[4]
        y(ksi, eta) = N[1](ksi, eta) * ye[1] + N[2](ksi, eta) * ye[2] + N[3](ksi, eta) * ye[3] + N[4](ksi, eta) * ye[4]

        # Creating jacobian Matrix for element e
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
        L = zeros(4, 4) # Local Mass Matrix
        
        # Precompute integrands for less allocations
        # Local potential Matrix V integrals
        integrands_V = [
            (ksi, eta) -> N[i](ksi, eta) * N[j](ksi, eta) * V_potential_func(x(ksi, eta), y(ksi, eta)) * abs(detJ(ksi, eta))
            for i in 1:4, j in 1:4
        ]

        # Compute local stiffness Matrix K integrals with αx and αy equal to 1 for the Shrodinger equation
        integrands_K = [
            (ksi, eta) -> gamma * (θNx[i](ksi, eta) * θNx[j](ksi, eta) + θNy[i](ksi, eta) * θNy[j](ksi, eta)) * abs(detJ(ksi, eta))
            for i in 1:4, j in 1:4
        ]

        # Compute local Mass Matrix integrals
        integrands_L = [
            (ksi, eta) -> N[i](ksi, eta) * N[j](ksi, eta) * abs(detJ(ksi, eta))
            for i in 1:4, j in 1:4
        ]

        # Add results to local Hamiltonian and Mass Matrix 
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
                va_H[index] = KV[i, j] # global value index at i, j of local Hamiltonian Matrix
                va_M[index] = L[i, j] # global value index at i, j of local Mass Matrix
                index = index + 1;
            end
        end
    end
    # Assembling global Hamiltonian Matrix = global stiff+ global potential, and Mass matrices using indices from loop
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
    println("Is H hermittian?", ishermitian(H))
    
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

    time1 = Base.time()

    time = time1 - time0
    time = round(time, digits = 2)
    println("Finished creating matrices from $noe elements in $time seconds")

    # Now our system looks like this: Aψ(n+1) = Βψn, almost ready for the solution with crank nicolson
    return A, B, lengthr, dt, boundary_nodes, tempo, step_size
end

# Creating inital wavefunction equation Ψ0 to solve with function solution
function wavefunction(x, y ; x0, y0, sigma, kx::Float64, ky::Float64)
    gaussian = exp(- ((x - x0)^2 + (y - y0)^2 ) / (2 * sigma^2))
    plane_wave = exp(-im * (kx * x + ky * y))
    return gaussian * plane_wave
end

# Creating function for getting solution from Ax=b where A is A,x is ψ(n+1) and b = B*ψn, also to test solvers etc
function solution(coords, nop, psi_zero, time, A, B, lengthr, dt, boundary_nodes, tempo, step_size, a, b, c, nx_half, ny_half, overlaps, iterations)
    
    # n steps of time
    time_domain = time[2] - time[1]
    n_steps = Int128(time_domain ÷ dt) + 1

    # Initialize psi_0 as a 1D vector
    psi_0 = Vector{ComplexF64}(undef, nop)

    # Extracting data from psi_zero function to psi_0 vector on our coordinate system
    for k in 1:nop
        psi_0[k] = psi_zero(coords[k,1], coords[k,2])
    end

    # Enforce Dirichlet boundary conditions (ψ = 0 at boundaries)
    psi_0[collect(boundary_nodes)] .= 0

    # Creating normalisation factor A
    A_ = sqrt(sum(psi_0' * psi_0 * step_size^2)) 

    psi_0 = (psi_0 / A_) # Normalising Ψ0
    temp = psi_0 # Keeping temp to plot inital function

    sum_check = sum(psi_0' * psi_0 * step_size^2) # Check if sum_check == 1 as it should be
    println("Normalized integral check: ", sum_check)

    # Calculate initial energy, (must stay the same until the end, since this is a closed system sim)
    E_initial = real(psi_0' * tempo * psi_0)
    println("Initial energy: $E_initial eV")


    # Initialising solution
    psi = similar(psi_0)
    A_LU = lu(A) # Preconditioning A for backslash

    # Time stepping loop

    # Keeping track of time for optimization
    t_start1 = Base.time()

    # Solving with \
    for n = 1:iterations

        # Solve psi
        psi = A_LU \ (B * psi_0)
        #println(sum(psi'*psi*step_size^2)) # Checking if integral stays the same
        # Asigning value to psi_0
        psi_0 = psi

    end
    temp1 = psi
    t_end1 = Base.time()

    # Reusing old solution to test my own solver
    psi_0 = temp
    
    println("Multigrid with domain decomposition method using $overlaps overlaps")
    t_start2 = Base.time()
    
    for n = 1:iterations

        # Solve psi
        psi = bicgstab_vic(A, B * psi_0, 1e-10, 150) # Bicgstab alone

        # Asigning value to psi_0
        psi_0 = psi

    end
    t_end2 = Base.time()

    # Reusing old solution to test solver
    psi_0 = temp
    
    t_start3 = Base.time()
    
    for n = 1:iterations

        # Solve psi
        psi = domain_decomposition(A, B * psi_0, a, b, c, nx_half, ny_half, overlaps, 1)

        # Asigning value to psi_0
        psi_0 = psi

    end
    t_end3 = Base.time()

    # Reusing old solution to test solver
    psi_0 = temp
    
    t_start4 = Base.time()
    
    for n = 1:iterations

        # Solve psi
        psi = domain_decomposition(A, B * psi_0, a, b, c, nx_half, ny_half, overlaps, 2)

        # Asigning value to psi_0
        psi_0 = psi

    end
    t_end4 = Base.time()

    # Final calculations

    final_E = real(psi_0' * tempo * psi_0)
    println("Energy: $final_E eV")

    time__ = t_end1 - t_start1
    time__ = round(time__, digits = 4)
    println("Time elapsed using backslash : $time__ seconds")
    
    time_2 = t_end2 - t_start2
    time_2 = round(time_2, digits = 4)
    println("Time elapsed using bicgstab solver alone: $time_2 seconds")

    time_3 = t_end3 - t_start3
    time_3 = round(time_3, digits = 4)
    println("Time elapsed using dd and bicgstab solver : $time_3 seconds")
    
    time4 = t_end4 - t_start4
    time4 = round(time4, digits = 4)
    println("Time elapsed using dd, multigrid and bicgstab solver : $time4 seconds")
    
    #println("Compare norms of psi solved with backslash $(norm(temp1)) vs multigrid $(norm(psi))")
    return psi, temp # return psi and the initial state 1 #
end

# Creating animated solution over desired time
function animated_solution(coords, nop, psi_zero, time, A, B, lengthr, dt, boundary_nodes, tempo, step_size, no_fr, no_afr)
    
    println("Iterating with Crank-Nicolson...")
    time0 = Base.time()
    # Initialize psi_0 as a 1D vector with a solution at each point so "nop"
    psi_0 = Vector{ComplexF64}(undef, nop)

    # Extracting data from psi_zero function to psi_0 vector on our coordinate system
    for k in 1:nop
        psi_0[k] = psi_zero(coords[k,1], coords[k,2])
    end

    # Enforce Dirichlet boundary conditions (ψ = 0 at boundaries)
    psi_0[collect(boundary_nodes)] .= 0

    # Creating normalisation factor A
    A_ = sqrt(sum(psi_0' * psi_0 * step_size^2)) 

    # Normalising Ψ0
    psi_0 = (psi_0 / A_)
    
    # Creating range for animation
    xg = yg = range(-1, 1, step = step_size)

    # Calculate global maximum for consistent scaling in plotting
    global_max = maximum(abs2.(psi_0)) * 1.8

    # Preconditioning
    A_LU = lu(A)

    # Initialisng frames and solution matrices to plot
    frames = []
    Z = fill(NaN, lengthr, lengthr)

    # Time evolution loop with Crank Nicolson
    for n in 1:no_fr # Produce

        # Solving psi with , right now \
        psi = A_LU \ (B * psi_0)

        # Time variable for plotting
        t = (n - 1) * dt * 1000
        t = round(t, digits = 2)
        # Capture every no_afrth frame
        if n % no_afr == 0 || n == 1 #800
            for k in 1:nop
                i = i = round(Int, (coords[k,:][2] + 1) / step_size) + 1
                j = round(Int, (coords[k,:][1] + 1) / step_size) + 1
                Z[i,j] = abs2(psi[k]) # ||^2 value of psi
            end
            # Plotting variable
            p = plot(
                surface(xg, yg, Z, # Surface for 3d plot
                    colormap = :viridis, # Matlab colormap viridis
                    colorbar = true, # Colorbar at the right to see at values
                    title = "Electron wavefunction |Ψ|² at t = $(round(t, digits=3)) femtoseconds",# Printing time as well
                    xlabel="x in nm", ylabel="y in nm", zlabel="|ψ(x, y)|²",
                    zlim = (0, global_max),
                    xlim = extrema(xg),# Extremes for x, y and the colorbar
                    ylim = extrema(yg),
                    clim = (0, global_max),
                    showscale = false,
                    camera = (30, 25, 1.0),
                    showlegend = false,
                    show_boundingbox = false,
                    overwrite_figure = false
                ),
                size = (800, 600),  # Resolution
                dpi = 300  # DPI
            )
            push!(frames, p)# Push p into frames for animation
        end
        psi_0 = psi # Ψ0 = Ψ to continue iterative method
    end
    time1 = Base.time()
    time = time1 - time0
    time = round(time, digits = 2)

    println("Finished iterating $no_fr frames in $time seconds")

    lengthf = length(frames)
    # Saving animation into variable
    println("Creating animation...")
    time0 = Base.time()
    anim = @animate for (i, frame) in enumerate(frames)
        plot(frame)  # Plot the current frame
        println("Progress: $i/$lengthf")  # Show current/total
    end

    # Printing energy level
    Energy = real(psi_0' * tempo * psi_0)
    println("Energy of closed system is: $Energy eV")
    time1 = Base.time()
    time = time1 - time0
    time = round(time, digits = 2)
    println("Animation created in $time seconds")
    return anim
end

# Creating 2 dimensional gauss_quadrature function to use to solve the integrals needed in FEM integrating from -1 to 1
function gauss_quad_2D(funct, n) 
    if n == 2
        ksi = (-0.5773502692, 0.5773502692)
        eta = (-0.5773502692, 0.5773502692)
        w = ones(2)
    elseif n == 3
        ksi = (0, -0.7745966692, 0.7745966692)
        eta = (0, -0.7745966692, 0.7745966692)
        w = (0.8888888889, 0.5555555556, 0.5555555556)
    elseif n == 4
        ksi = (-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116)
        eta = (-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116)
        w = (0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451)
    else
        error("Only n = 2, 3, 4 are supported")
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
function V_function(V_flag::Int64, V0 = 0.0, x0 = 0.0, y0 = -0.5, r = 0.2)
    if V_flag == 1 # Infinite potential well with bottom at 0.0 or at selected V0
        return (x,y) -> V0 # Which equals 0.0 at default
    elseif V_flag == 2 # Barrier of potential V0 at area around x0 and y0"default center at (-0.3, -0.5) with radius 0.1"
        return (x, y) -> ((x - x0)^2 + (y - y0)^2 <= r^2 ? V0 : 0.0)
    elseif V_flag == 3 # Well of potential V0
        return (x, y) -> ((x - x0)^2 + (y - y0)^2 <= r^2 ? -V0 : 0.0)
    end
end

# Solver function Bicgstab since system has non spd and complex matrices A and B * psi_0
function bicgstab_vic(A, b, tol, Nmax)

    # Getting size of A
    n = size(A, 1)

    # Preallocate vectors
    x = zeros(ComplexF64, n)
    r = similar(b) # Residual vector
    r .= b .- A * x
    r_hat = copy(r)
    p = copy(r) # Vector p
    v = similar(b)
    s = similar(b)
    s_hat = similar(b)
    t_hat = similar(b)
    t = similar(b) 
    ro = dot(r_hat, r) # ρ0

    # Iteration counter for optimization
    iter = 0

    # b Norm
    norm_b = norm(b)

    # Using preconditioner, as a test using ilu
    K = ilu(A, τ = 1e-6)

    for i in 1:Nmax
        mul!(v, A, p)
        alpha = ro / dot(r_hat, v)
        h = x + alpha * p
        s = r - alpha * v
        s_res = norm(s) / norm_b

        # Break if reached tolerance
        if s_res < tol
            x = h
            #println("Converged at iteration $iter with residual norm: $s_res")
            break
        end
        s_hat .= K \ s
        t = A * s_hat
        t_hat .= K \ t
        omega = dot(t_hat, s_hat) / dot(t_hat, t_hat)
        x = h + omega * s_hat
        r = s - omega * t

        rel_res = norm(r) / norm_b
        
        if rel_res < tol
            #println("Converged at iteration $iter with residual norm: $rel_res")
            break
        end
        temp = ro
        ro = dot(r_hat, r)
        beta = (ro / temp) * (alpha / omega)
        p = r + beta * (p - omega * v)
        iter += 1
        if i == Nmax
            #println("Exiting with max iterations..")
        end
    end
    return x
end

# Function for domain decomposition method
function domain_decomposition(A, b, s1, s2, s3, nx_half, ny_half, overlap, flag) # Matrix A and sub matrices from grid A, b, 
    
    # error handling
    if overlap > nx_half - 1 || overlap > ny_half - 1 # n_half = nx / 2 + 1, thats why -1 is there
        println("Overlap needs to  less than or equal to n / 2 ...")
        return 0
    end

    # Getting size of A
    n = size(A, 1)

    # Preallocate solution and residual vectors
    x = zeros(ComplexF64, n, 1)
    r = similar(b)

    # s1, s2, s3 modification to get inner and boundary nodes, inner2 = s2
    inner1 = s1[1: end - 1, :] # Last line is the boundary line same as first line of s2
    inner3 = s3[:, 2:end] # First collumn is the boundary collumn aligning with s2

    # Overlap matrices
    over1 = s2[1:overlap, :]
    over21 = inner1[end - overlap + 1:end, :]
    over23 = inner3[:, 1:overlap]
    over3 = s2[:, end - overlap + 1:end]

    # Indices for submatrices
    inds1 = vec(transpose(vcat(inner1, over1)))
    inds2 = vec(hcat(transpose(vcat(over21, s2)), over23))
    inds2 = sort(inds2)
    inds3 = vec(transpose(hcat(over3, inner3)))

    # Sumbatrices of A
    A1 = A[inds1, inds1]
    A2 = A[inds2, inds2]
    A3 = A[inds3, inds3]
    # Resulting matrices will be of size nop x nop where nop is the number of nodes within their subdomain, they will be square

    # Creating index Matrix that gets wanted solution from x2 after solving, which is where s2 lies, meaning, we dont want the solution
    # that is on over21 and over23 but on top of s2
    d = zeros(Int, length(vec(s2)), 1)
    index = length(vec(over21)) + 1 # index starts from where over21 ends, also used to give value to d
    indexx = 1
    for k = 1:size(s2)[1] # Row loop
        for i = 1: size(s2)[2]
            d[indexx] = index
            index += 1
            indexx +=1
        end
        index +=  size(over23)[2] # skip over23
    end

    # Same for s3 to not include the first collumn of the solution
    c = zeros(Int, length(vec(inner3)), 1)
    index = overlap + 1 # skips first collumn
    indexx = 1 
    for k = 1:size(inner3)[1] # Row loop
        for i = 1: size(inner3)[2]
            c[indexx] = index
            index += 1
            indexx +=1
        end
        index += overlap # skips first collumn
    end
    
    # Max iterations and tolerance
    Nmax = 200;
    tolDD = 10e-10
    tolMG = 10e-4
    tolBic = 10e-5
    it = 20

    if flag == 1 # Solving with bicgstab
        for i = 1:Nmax
            # Residual
            r = b - A * x;

            # Convergence criterion and fail print
            if(norm(r)/norm(b)<tolDD)
                #println("Converged in $i iterations with residual ", norm(r)/norm(b))
                break
            end
            if i == Nmax
                println("Reached maximum iterations without convergence")
            end

            # Solve all systems with residual in selected indices
            
            # Solve first system
            
            x1 = bicgstab_vic(A1, r[inds1], tolBic, it) # multigrid_vic(nx_1, ny_1, overlap, A1, r[inds1], 3, 2, tol, 1) #
            x1 = x1[1:(length(inner1))] # Stays the same since indexing is the same

            # Solve second system
            x2 = bicgstab_vic(A2, r[inds2], tolBic, it) # multigrid_vic(nx_2, ny_2, overlap, A2, r[inds2], 3, 2, tol, 2) #
            x2 = x2[d]

            # Solve the third system
            x3 = bicgstab_vic(A3, r[inds3], tolBic, it) # multigrid_vic(nx_3, ny_3, overlap, A3, r[inds3], 3, 2, tol, 1) #
            x3 = x3[c]
    
            # Update solution
            x[sort!(vec(inner1))] += x1
            x[sort!(vec(s2))] += x2
            x[sort!(vec(inner3))] += x3
        end
    elseif flag == 2 # Using multigrid and bicgstab as a smoother

        # Parameters
        n1 = 6
        n2 = 4

        # Getting sizes for multigrid

        nx_1 = size(inner1)[2]
        ny_1 = size(inner1)[1] + overlap

        nx_2 = size(s2)[1] + overlap
        ny_2 = nx_2

        nx_3 = size(inner3)[2] + overlap
        ny_3 = size(inner3)[1]

        time_start = Base.time()
        # Preallocating matrices for multigrid
        A1_, R1, P1, levels1 = multigrid_vic_hierarchy(nx_1, ny_1, overlap, A1)
        A2_, R2, P2, levels2 = multigrid_vic_hierarchy(nx_2, ny_2, overlap, A2, 2)
        A3_, R3, P3, levels3 = multigrid_vic_hierarchy(nx_3, ny_3, overlap, A3)

        time_end = Base.time()
        time = time_end - time_start
        time = round(time, digits = 4)
        println("Time elapsed creating matrices for multigrid before loop:", time)

        for i = 1:Nmax
            # Residual
            r = b - A * x;

            # Convergence criterion and fail print
            if(norm(r)/norm(b)<tolDD)
                #println("Converged in $i iterations with residual ", norm(r)/norm(b))
                break
            end
            if i == Nmax
                println("Reached maximum iterations without convergence")
            end

            # Solve all systems with residual in selected indices
            
            # Solve first system
            
            x1 = multigrid_vic(A1_, R1, P1, levels1, r[inds1], n1, n2, tolMG)
            x1 = x1[1:(length(inner1))] # Stays the same since indexing is the same

            # Solve second system
            x2 = multigrid_vic(A2_, R2, P2, levels2, r[inds2], n1, n2, tolMG)
            x2 = x2[d]

            # Solve the third system
            x3 = multigrid_vic(A3_, R3, P3, levels3, r[inds3], n1, n2, tolMG)
            x3 = x3[c]
            
            # Update solution
            x[sort!(vec(inner1))] += x1
            x[sort!(vec(s2))] += x2
            x[sort!(vec(inner3))] += x3
        end
    end
    return x
end
function multigrid_vic_hierarchy(nox, noy, overlap, A1, flag = 1)
    
    #  First I will be creating the restriction and prolongation matrices by using the weighted average of nine fine grid values to create
    # for non_f = fine grid number of nodes in a direction, to find the coarse nodes non_c = floor((non_f - 1) / 2) + 1 so that if the interval number is odd it becomes even
    # example 5x x 6y nodes -> 4x x 5y intervals, non_cx = 3 and non_cy = 3 also, containing the information of all the previous finer nodes


    # Getting number of nodes for each coarser level

    levels = 1
    nx = zeros(Int, 21, 1) # A size of 21 is sufficient to have as many as 2 million nodes, set number less expensive than calculating the levels
    ny = similar(nx)

    nx[1] = nox # Fine grid values at level 1
    ny[1] = noy

    #println("x: $(nx[1]), y: $(ny[1]) at level: $levels")

    for i = 2: 21
        nx[i] = floor((nx[i - 1] - 1) / 2) + 1
        ny[i] = floor((ny[i - 1] - 1) / 2) + 1
        levels += 1
        #println("x: $(nx[i]), y: $(ny[i]) at level: $levels")
        if nx[i] <= 3 || ny[i] <= 3
            break
        end
    end
    levelss = 1
    if flag == 2
        n_inner = zeros(Int, 21, 1) 
        n_inner[1] = nox - overlap
        #println("Inner nodes: $(n_inner[1]) at level: $levelss")
        for i = 2: 21
            n_inner[i] = floor((n_inner[i - 1] - 1) / 2) + 1
            levelss += 1
            #println("Inner nodes: $(n_inner[i]) at level: $levels")
            if nx[i] <= 3
                break
            elseif n_inner[i] <= 3
                for j = i + 1: levels
                    n_inner[j] = n_inner[i]
                    levelss += 1
                    #println("Inner nodes: $(n_inner[j]) at level: $levelss")
                end
                break
            end
        end
    end
    # Create restriction and prolongation matrices which are for levels 2-total number of levels since the first level doesnt need them
    
    R = Vector{SparseMatrixCSC{Float64, Int}}(undef, levels - 1) # Vectors of sparse matrices
    P = Vector{SparseMatrixCSC{Float64, Int}}(undef, levels - 1)

    for i = 2:levels
        if flag == 1 # For A1 and A3 which are rectangles
            R[i - 1] = build_Restriction_matrix(flag, nx[i - 1], nx[i], ny[i - 1], ny[i])
        else # If its A2 and it has an L shape
            R[i - 1] = build_Restriction_matrix(flag, nx[i - 1], nx[i], ny[i - 1], ny[i], n_inner[i - 1], n_inner[i])
        end
        # Create prolongation matrix from the restriction
        P[i - 1] = 4 * R[i - 1]'
    end 
    # Inniate A

    A = Vector{SparseMatrixCSC{ComplexF64, Int}}(undef, levels) # Vector of sparse matrices

    # Build A's top level (finest)
    A[1] = A1

    # Use Galerkin method to build A on each level
    for i = 2: levels
        A[i] = R[i - 1] * A[i - 1] * P[i - 1]
    end

    return A, R, P, levels
end
# Creating multigrid function with V cycles for my 2D grid, as it the most efficient, to use as preconditioner along with domain decomposition
#so this multigrid will be tailored for rectangle grids (not L shape since we use decomposition)
function multigrid_vic(A, R, P, levels, B, n1, n2, tol)

    # Parameters
    Nmax = 17

    # Iniate b and x vectors 
    b = Vector{Vector{ComplexF64}}(undef, levels)
    x = similar(b)

    for i = 1: levels
        b[i] = zeros(ComplexF64, size(A[i])[1])
        x[i] = similar(b[i])
    end
    # Build the RHS's top level
    b[1] = B

    # Norm of RHS
    norm_b = norm(b[1])

    # Iterate over all levels
    for it = 1: Nmax
        for i = 2: levels - 1
            x[i] = bicgstab_vic(A[i], b[i], 10e-10, n1)
            b[i + 1] = R[i] * (b[i] - A[i] * x[i])
            x[i + 1][:].= 0
        end
        x[end] = A[end] \ b[end] 

        for i = levels - 1 : -1: 1
            x[i] = x[i] + P[i] * x[i+1]
            x[i] = bicgstab_vic(A[i], b[i], 10e-10, n2)
        end

        nrm = norm(b[1] - A[1] * x[1])
        if nrm < norm_b * tol
            break;
        end
        if it == Nmax
            println("Exiting with max iterations on multigrid")
            break;
        end
    end
    return x[1]
end

# Map function to get active node mappings for L-shape and not get the top right corner.
function map(flag, nx::Int, n_inner::Int = 0, ny::Int = nx)

    map_index = zeros(Int, ny, nx) # Saves the indicies of the active nodes
    map_coords = Tuple{Int, Int}[] # Stores the coordinates of active node index

    n_act = 0 # Active number of nodes

    # Loop to check if the node is active and get the indices and coordinate
    for y = 1:ny
        for x = 1:nx
            if ((x <= n_inner && y <= ny) || (y <= n_inner && x <= nx)) && flag == 2
                n_act += 1
                map_index[y, x] = n_act
                push!(map_coords, (y, x))
            elseif flag == 1
                n_act += 1
                map_index[y, x] = n_act
                push!(map_coords, (y, x))
            end
        end
    end
    return n_act, map_index, map_coords
end

# Function to build the L-shaped Restriction Matrix using fine and coarse grid bounding box and inner box parameters
function build_Restriction_matrix(flag, nfx_bb::Int, ncx_bb::Int, nfy_bb::Int = nfx_bb, ncy_bb::Int = ncx_bb, nf_inner::Int = 0, nc_inner::Int = 0)
    if flag == 2
        nf_act, f_map_index, _ = map(flag, nfx_bb, nf_inner) # Getting the fine active nodes as well as the fine map indices
        nc_act, _, c_map_coords = map(flag, ncx_bb, nc_inner) # Getting the coarse active nodes as well as the coordinates for the coarse map
    else
        nf_act, f_map_index, _ = map(flag, nfx_bb, 0, nfy_bb) # Getting the fine active nodes as well as the fine map indices
        nc_act, _, c_map_coords = map(flag, ncx_bb, 0, ncy_bb) # Getting the coarse active nodes as well as the coordinates for the coarse map
    end
    if nc_act == 0 || nf_act == 0
        return spzeros(Float64, nc_act, nf_act)
    end

    I = Int[] # Row indices (active coarse node L-indices)
    J = Int[] # Collumn indices (active fine node L-indices)
    V = Float64[] # Stencil weights

    stencil = 1/16 * [1 2 1; 2 4 2; 1 2 1] # 9 node stencil, diagonal 1/16, horizontal and vertical 1/8 and central 1/4 adding up to 1

    for  kcL = 1: nc_act # Iterate over each active coarse node
        yc_bb, xc_bb = c_map_coords[kcL] # 2D bounding box coords of this coarse node on y, x

        # Getting central fine grid point correspoing to where the coarse node is
        yf_ctr = 2 * yc_bb - 1 # Fine grid center y 
        xf_ctr = 2 * xc_bb - 1 # Fine grid center x
        
        for s_y = 1:3 # Stencil in the y-dim (1,2,3)
            for s_x = 1:3 # Stencil in x-dim
                
                weight = stencil[s_y, s_x]
                if weight == 0.0 continue end # If there is no node continue to the next node

                off_y = s_y - 2 # Stencil offset (-1,0,1)
                off_x = s_x - 2 # Stencil offset (-1,0,1)
                
                tfy = yf_ctr + off_y # Target fine y
                tfx = xf_ctr + off_x # Target fine x

                if (1 <= tfy <= nfy_bb) && (1 <= tfx <= nfx_bb)
                    kfL = f_map_index[tfy, tfx]
                    if kfL>0
                        push!(I, kcL)
                        push!(J, kfL)
                        push!(V, weight)
                    end
                end
            end
        end
    end
    
    if isempty(I)
        return spzeros(Float64, nc_act, nf_act)
    end

    R = sparse(I, J, V, nc_act, nf_act)
    #=row_sums = sum(R, dims=2)
    for i in 1:nc_act
        if row_sums[i] != 0
            R[i, :] = R[i, :] ./ row_sums[i]
        end
    end=#
    return R
end
end # Module end