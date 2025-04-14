#= Victor Emmanuel Meimaris 23/3/25: Creating Module with all functions used on main.jl, I will be
explaining the numbers in my report. =#
module Functions
export grid, righthandside, gauss_quad_2D
using SparseArrays
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
            l2g[index, :] = [a[i, j], a[i + 1, j], a[i + 1, j + 1], a[i, j + 1]]
            index +=1
        end
    end
    # For bottom rectangle
    for i = 1:(num_of_squares_y ÷ 2)
        @inbounds for j = 1: num_of_squares_x
            l2g[index, :] = [b[i, j], b[i + 1, j], b[i + 1, j + 1], b[i, j + 1]]
            index += 1
        end
    end
    return coords, l2g ,noe, nop, boundary_nodes, max_length
end
# Creating function which will be returning the stiffness matrix and the right hand side of the system with F.E.M
function matrix_equation(funct_g, coords::Matrix{Float64}, l2g::Matrix{Int64}, noe::Int64, nop::Int64, boundary_nodes::Set, max_length, V_flag::Int64) # Can just put one variable results = Functions.grid(domain, max_length) in main
    
    # For this script I will assume Neuman and robin coefficients to be zero and then might integrate them for certain problems like scattering
    # I will also assume Dirichlet boundary conditions where the probability ψ is zero in the boundaries, as in this general function I want to model
    # a closed system. If an open system would be simulated, the matrix p would be needed as well, whichs contains the neuman and robin coefficients

    # Assemble stiffness matrix indices (triplets for sparse matrix)
    ia = zeros(noe * 16); # Row index, noe * 16 because the two for loops for each square = 16 iterations for each element --> noe * 16
    ja = zeros(noe * 16); # Column index
    va = zeros(noe * 16); # Value index

    # Initialize global right hand side vector b
    b = zeros(nop, 1);

    # Shape/basis functions for integration
    N = [
        (ksi,eta) -> 1/4 * (1 - ksi) * (1 + eta),# N1
        (ksi,eta) -> 1/4 * (1 - ksi) * (1 - eta),# N2
        (ksi,eta) -> 1/4 * (1 + ksi) * (1 - eta),# N3
        (ksi,eta) -> 1/4 * (1 + ksi) * (1 + eta)# N4
    ]

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

        # Initialise local right hand side vector f and assemble global right hand side vector b
        fl = 0
        # Precompute integrands for less allocations
        integrands_fl = [
            (ksi, eta) -> N[i](ksi, eta) * funct_g(x(ksi, eta), y(ksi, eta)) * detJ(ksi, eta)
            for i in 1:4
        ]
        @inbounds for i = 1:4
            fl = gauss_quad_2D(integrands_fl[i], 2)
            # Assemble global right hand side vector b
            b[l2g[e, i]] += fl
        end

        # Define θNx and θNy as arrays of functions, will not include 1/detJ(ksi, eta) since it is multiplied by the same function anyway so i will just / in the integrand 
        #for efficiency instead of calculating it 8 times here and / it another 48 times
        θNx = [
            (ksi, eta) -> 0.25 * (-J22(ksi) * (eta + 1) - J12(eta) * (1 - ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (J22(ksi) * (eta - 1) + J12(eta) * (1 - ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (-J22(ksi) * (eta - 1) + J12(eta) * (ksi + 1)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (J22(ksi) * (eta + 1) - J12(eta) * (ksi + 1)) / detJ(ksi, eta)
        ]
        θNy = [
            (ksi, eta) -> 0.25 * (J21(ksi) * (1 + eta) + J11(eta) * (1 - ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (J21(ksi) * (1 - eta) - J11(eta) * (1 - ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (-J21(ksi) * (1 - eta) - J11(eta) * (1 + ksi)) / detJ(ksi, eta),
            (ksi, eta) -> 0.25 * (-J21(ksi) * (1 + eta) + J11(eta) * (1 + ksi)) / detJ(ksi, eta)
        ]
       
        # Local stiffness matrix M, mass matrix not needed
        M = zeros(4, 4)
        # Compute local stiffness matrix with αx and αy equal to 1 for this project, for more advanced application this could change
        for i = 1:4
            @inbounds for j = 1:4
                integrand_M(ksi, eta) = -(θNx[i](ksi, eta) * θNx[j](ksi, eta) + θNy[i](ksi, eta) * θNy[j](ksi, eta)) * detJ(ksi, eta)
                M[i, j] = gauss_quad_2D(integrand_M, 3) # value of local stiffness matrix at i,j
            end
        end

        # Assemble global stiffness matrix indices from the local element matrices 
        for i = 1:4
            @inbounds for j = 1:4
                ia[index] = l2g[e, i] # global node index at row i of local matrix
                ja[index] = l2g[e, j] # global node index at collumn j
                va[index] = M[i, j] # global value index at i, j of local stiffness matrix
                index = index + 1;
            end
        end
    end
    # Assembling global stiffness matrix using indices
    K = sparse(ia, ja, va)

    # Enforce dirichlet boundary conditions
    boundary_nodes = collect(boundary_nodes)  # Convert Set to Vector
    K[boundary_nodes, :] .= 0
    b[boundary_nodes, :] .= 0
    for node in boundary_nodes
        K[node, node] = 1
    end
    return K, b
end
# Creating function for getting solution from
function solution(K, b)
    q = K\b
    return q
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
function V_potential(V_flag::Int64)
    q = 0 # Flux term for neumman 
    gamma = 0 # Robin term for mixed conditions
    if V_flag == 1 # Infinite potential well
        return 1
    elseif V_flag == 2
        q = 1
        gamma = 1
    end
    return q, gamma
end
end # Module end