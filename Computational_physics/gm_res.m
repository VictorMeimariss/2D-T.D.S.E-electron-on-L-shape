function [x, iter, flag] = gm_res(A, b, m, max_restarts, tol)
    flag = 0;
    % Getting size of matrix A 
    n = size(A, 1);
    iter = 0;
    x = zeros(n, 1);
    for restart = 1:max_restarts

        r = b - A * x; % Compute residual each time gmres restarts

        [V, H, m] = arnoldi(A, r, m); % Use arnoldi method to find V and H

        beta = norm(r);
        e1 = zeros(m+1, 1);
        e1(1, 1) = beta; % This is β * e1
        y = H(1:m, :) \ e1(1:m); %y = H \ e1; % Solve y

        x = x + V(:, 1:m) * y; % Vm = V nxm we basically delete the last collumn
        r = b - A * x;
        iter = iter + m;
        if norm(r) < tol
            flag = 1;
            return;
        end
    end
end
% Creating arnoldi function to get basis matrix with u vectors and Hessenberg matrix
function[V, H, m] = arnoldi(A, r, m)
    
    % Getting size of matrix A 
    n = size(A,1);
    % Allocating memory to H(m+1 x m)
    H = zeros(m+1, m);
    % Allocating memory to V(n x m+1)
    V = zeros(n, m+1);

    V(:, 1) = r/ norm(r); % The first column of V consists of the normalised ||u||2=1 vectors
    for j = 1:m
        w = A * V(:, j);
        for i = 1:j
            H(i, j) = V(:, i)' * w;
            w = w - H(i, j) * V(:, i);
        end
        H(j+1, j) = norm(w);
        if abs(H(j+1, j)) < sqrt(eps)
            H = H(1:j, 1:j);
            V = V(:, 1:j);
            m = j;
            return;
        end
        V(:, j+1) = w / H(j+1, j);
    end
end