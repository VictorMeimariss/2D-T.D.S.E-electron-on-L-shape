include("functions.jl")  # Include the file
import .Functions  # Use the module (note the dot '.')

domain = (-1, 1)
max_length = 1.0

results = Functions.grid(domain, max_length)  # Get tuple 
results_1 = Functions.righthandside(results...) # Splat the tuple
println(results_1) # Have to create seperate function for solution and seperate for grid
#f = (x, y) -> x^2 + y^2  # Anonymous function
#approx = Functions.gauss_quad_2D(f, 2)
#println("Approximation for n=2: ", approx)

#approx = Functions.gauss_quad_2D((θ, η) -> sin(θ) + cos(η), 3)  # Inline anonymous function
#println("Approximation for n=3: ", approx)