include("functions.jl")  # Include the file
import .Functions  # Use the module (note the dot '.')

domain = (-1, 1)
max_length = 1.0

results = Functions.grid(domain, max_length)  # Call the function
println(results) # Have to create seperate function for solution and seperate for grid