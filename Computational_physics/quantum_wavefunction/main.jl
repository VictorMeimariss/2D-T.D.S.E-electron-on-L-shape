using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Plots

# 1️⃣ Define the domain (L-shaped region)
domain = (0, 1, 0, 1)  # A square domain for now (we'll make L-shape later)

# 2️⃣ Generate a rectangular mesh
model = CartesianDiscreteModel(domain, (20, 20))  # 20x20 grid

# 3️⃣ Define function space
order = 1  # Linear basis functions (FEM order)
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe, conformity=:H1)

# 4️⃣ Define the PDE: Poisson equation
Ω = Triangulation(model)
dΩ = Measure(Ω, 2 * order)

# Weak form of the Poisson equation: ∫(∇u ⋅ ∇v) dΩ = ∫(f ⋅ v) dΩ
a(u, v) = ∫(∇(u) ⋅ ∇(v)) * dΩ  # Bilinear form
b(v) = ∫(1 * v) * dΩ  # Linear form (source term f = 1)

# 5️⃣ Solve the FEM system
op = AffineFEOperator(a, b, V, V)
uh = solve(op)

# 6️⃣ Extract nodal values for plotting
# Use Gridap's `get_cell_coordinates` to extract nodal coordinates
cell_coords = get_cell_coordinates(Ω)
x = vcat(cell_coords...)  # Flatten the list of coordinates

# Convert VectorValue{2, Float64} to separate arrays for x and y coordinates
x_coords = [point[1] for point in x]
y_coords = [point[2] for point in x]

# Evaluate the solution at the nodes
uh_values = uh.(x)

# 7️⃣ Plot the solution
plotlyjs()
scatter(x_coords, y_coords, uh_values, marker_z=uh_values, label="Solution", color=:viridis)
xlabel!("x")
ylabel!("y")
title!("FEM Solution for Poisson Equation")