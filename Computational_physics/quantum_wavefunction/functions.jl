# Victor Emmanuel Meimaris 23/3/25: Creating Module with all functions used on main.jl

module Functions
export coords
# Creating Coordinate function
function coords(domain::Tuple{Real, Real}, max_length::Float64) # num_of_squares_x::UInt16, num_of_squares_y::UInt16

    # Define domain length in any axis since x and y have the same domain
    domain_length = domain[2] - domain[1]

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
    accuracy_level = 2 # Might change these later, put it as input or remove it
    num_of_squares_y = temp * accuracy_level
    num_of_squares_x = temp * accuracy_level # This could be different to y if multiple of 2

    # Create square parts of the global matrix and then merge them with v and h cat into one global matrix g
    a = reshape(1:((num_of_squares_x ÷ 2) + 1) * ((num_of_squares_y ÷ 2) + 1), (num_of_squares_x ÷ 2) + 1, (num_of_squares_y ÷ 2) + 1)
    a = transpose(a)

    # Creating (n+1, n) sparse matrix for domain outside of Omega and boundaries with it to merge with a
    b = zeros(Int, (num_of_squares_y ÷ 2) + 1, num_of_squares_x ÷ 2)
    temp = a[end, end] + 1
    for i = 1:num_of_squares_x ÷ 2
        b[end, i] = temp
        temp +=1
    end

    # Merging a and b horizontally to form half of global matrix
    g = hcat(a, b) 

    # Number of elements and nodes
    noe = ((num_of_squares_y ÷ 2) * num_of_squares_x) + ((num_of_squares_y ÷ 2) * (num_of_squares_x ÷ 2))  # Number of elements
    nop = (num_of_squares_x + 1) * (num_of_squares_y + 1) - length(b) + num_of_squares_x ÷ 2 # Number of nodes
    
    # Creating rest of the global map and merging it with g
    # c is the remainder and goes from temp to nop - length(b) + num_of_squares_y ÷ 2 
    c = reshape(temp:nop, num_of_squares_x + 1, num_of_squares_y ÷ 2)
    c = transpose(c)

    # Merge c and g into g the global matrix
    g = vcat(g, c)

    # Generate coordinates
    x = range(domain[1], domain[2], num_of_squares_x + 1)
    y = range(domain[1], domain[2], num_of_squares_y + 1)

    # Define local to global map " l2g "
    l2g = zeros(Int, noe, 4) # Each element has 4 nodes
    temp = 1
    # For top left corner of domain
    for i = 1: num_of_squares_y ÷ 2
        for j = 1: num_of_squares_x ÷ 2
            l2g[temp, :] = [g[i, j], g[i + 1, j], g[i + 1, j + 1], g[i, j + 1]]
            temp += 1
        end
    end
    # For bottom half
    for i = (num_of_squares_y ÷ 2) + 1: num_of_squares_y
        for j = 1: num_of_squares_x
            l2g[temp, :] = [g[i, j], g[i + 1, j], g[i + 1, j + 1], g[i, j + 1]]
            temp += 1
        end
    end
    return g, l2g
end
end # Module end