clc;
clear;

% Domain size
L = [0, 1];

% Define anonymous force function of PDE on the right-hand side
f = @(x, y) 2 * pi^2 * sin(pi * x) .* sin(pi * y); % Use element-wise multiplication

% Define number of squares in each direction
num_of_squares_x = 2;
num_of_squares_y = 2;

% Global Matrix from 1-9 for ex: 123; 456; 789 indicating global nodes
a = reshape(1:(num_of_squares_x + 1) * (num_of_squares_y + 1), num_of_squares_x + 1, num_of_squares_y + 1).';

% Generate x and y coordinates 
x = linspace(L(1), L(2), num_of_squares_x + 1);
y = linspace(L(1), L(2), num_of_squares_y + 1);

% Number of elements and nodes
noe = 2 * num_of_squares_x * num_of_squares_y; % Number of elements (triangles)
nop = (num_of_squares_x + 1) * (num_of_squares_y + 1); % Number of nodes

% Define local to global map
l2g = zeros(noe, 3); % Each element has three nodes
index = 1;
for i = 1:num_of_squares_y
    for j = 1:num_of_squares_x
        % First triangle in square
        l2g(index, :) = [a(i, j), a(i, j + 1), a(i + 1, j)];
        % Second triangle in square
        l2g(index + 1, :) = [a(i + 1, j + 1), a(i + 1, j), a(i, j + 1)];
        index = index + 2;
    end
end

% Coordinates
[xc, yc] = ndgrid(x, y); % Creates two matrices: one with each line as x coordinates ("xc") and one as each column as y coordinates ("yc")
coords = [xc(:), yc(:)]; % Two-column matrix with the first column having x coordinates and the second having y coordinates

% Assemble stiffness matrix indices (triplets for sparse matrix)
ia = zeros(noe, 1); % Row index
ja = zeros(noe, 1); % Column index
va = zeros(noe, 1); % Value index

% Initialize global force vector
F = zeros(nop, 1);

% Shape(basis ) functions for integration
N1 = @(ksi, hta) 1 - ksi - hta;
N2 = @(ksi, hta) ksi;
N3 = @(ksi, hta) hta;

% Iterate over elements to find global stiffness and force matrices
index = 1;
for e = 1:noe
    xe = coords(l2g(e, :), 1);
    ye = coords(l2g(e, :), 2);

    % Compute edge vectors
    y23 = ye(2) - ye(3);
    x32 = xe(3) - xe(2);
    y31 = ye(3) - ye(1);
    x13 = xe(1) - xe(3);
    x21 = xe(2) - xe(1);
    y21 = ye(2) - ye(1);

    x31 = - x13;

    % Compute local stiffness matrix
    Ae = (x21 * y31 - x13 * y21) / 2; % Triangle area that's why /2
    M = (1 / (4 * Ae)) * ...
        [y23^2 + x32^2, y23 * y31 + x32 * x13, y23 * y21 + x32 * x21;
         y23 * y31 + x32 * x13, y31^2 + x13^2, y31 * y21 + x13 * x21;
         y23 * y21 + x32 * x21, y31 * y21 + x13 * x21, y21^2 + x21^2];

    % Assemble global stiffness matrix from the local element matrices
    for i = 1:3
        for j = 1:3
            ia(index) = l2g(e, i); % global node index at row i of local matrix
            ja(index) = l2g(e, j); % global node index at collumn j
            va(index) = M(i, j); % value of local stiffness matrix at i,j
            index = index + 1;
        end
    end

    % Compute Jacobian determinant to transform integrals from the global to local coordinates
    J = [x21, x31;
         y21, y31];
    detJ = abs(det(J));

    % Define mapping from local (ksi, hta) to global (x, y) coordinates
    % This done by multiplying with basis functions
    %( I am using function overloading, don't confuse x,y with the previous x, y)
    x = @(ksi, hta) xe(1) * N1(ksi, hta) + xe(2) * N2(ksi, hta) + xe(3) * N3(ksi, hta);
    y = @(ksi, hta) ye(1) * N1(ksi, hta) + ye(2) * N2(ksi, hta) + ye(3) * N3(ksi, hta);

    % Compute local force vector
    Fl = zeros(3, 1);
    for i = 1:3
        Ni = {N1, N2, N3};
        Fl(i) = integral2(@(ksi, hta) f(x(ksi, hta), y(ksi, hta)) .* Ni{i}(ksi, hta) .* detJ, 0, 1, 0, @(ksi) 1 - ksi);
    end

    % Assemble global force vector
    for i = 1:3
        F(l2g(e, i)) = F(l2g(e, i)) + Fl(i);
    end
end

% Assemble global stiffness matrix
K = sparse(ia, ja, va);

% Identify boundary nodes
boundary_nodes = unique([a(1, :), a(end, :), a(:, 1)', a(:, end)']); % All boundary nodes

% unique ensures that duplicate nodes are only counted once ( corner nodes) 

% To enforce boundary conditions u = 0 on (x, y) ∈ dΩ which means 0 on the
% boundaries, we firstly set K to zero on its boundaries ( which would mean
% that for example in a 9x9 the central node would be unaffected), and then
% to avoid singularity when solving the linear equation of boundary nodes
% turn the diagonal node's value to 1(except for the nodes that aren't affected
% by the boundary conditions). The reason we even zero out nodes 


% Enforce Dirichlet boundary conditions
K(boundary_nodes, :) = 0; % This ensures that the equation for the boundary node
% i does not depend on the values of other nodes
K(sub2ind(size(K), boundary_nodes, boundary_nodes)) = 1; % Setting K=1 diagonally on boundary nodes to avoid singularity
F(boundary_nodes) = 0; % F = 0 on the boundaries


% Solve the linear system
q = K \ F;

% Display the solution
%disp('Solution q:');
%disp(q);


% Visualize solution on triangular mesh with trisurf
figure;
trisurf(l2g, coords(:, 1), coords(:, 2), q);
title('FEM Solution on Triangular Mesh');
xlabel('x');
ylabel('y');
zlabel('q(x, y)');
colorbar; % To understand what value each color represents

% Display error, the analytical solution is u(x,y) = sin(pi * x) .* sin(pi * y) 
g = @(x,y) sin(pi * x) .* sin(pi * y);
x = linspace(L(1), L(2), num_of_squares_x + 1);
y = linspace(L(1), L(2), num_of_squares_y + 1);

error = zeros(nop,1);
for i=1:num_of_squares_x + 1
    for j=1:num_of_squares_y + 1
        node_index = (num_of_squares_y + 1 - j) * (num_of_squares_x + 1) + (num_of_squares_x + 1 - i) + 1;
        analytical_value = g(x(i), y(j));
        error(node_index) = abs(analytical_value - q(node_index));
    end
end
figure;
plot(error);
xlabel('Node index');
ylabel('Absolute error');