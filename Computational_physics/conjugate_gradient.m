% Conjugate gradient function
function [x, iter, flag] = conjugate_gradient(A, b, tol)
    
    % Getting size of matrix A
    n = size(A, 1);

    % Solution
    x = zeros(n, 1);

    % Residual
    r = b - A * x;
    p = r;

    flag = 0;

    % Norm of rbs
    norm0 = norm(r);
    
    % Max iterations
    Nmax = 1000; 

    % Iteration counter
    iter = 0;

    for i = 1 : Nmax 
        
        rr = r' * r; 
        r0 = r;
        Ap = A * p;
        a = rr / (p' * Ap);
    
        x = x + a * p; % x(n+1) = x(n) + a * p 
        r = r - a * Ap;
    
        if norm(r) < tol * norm0
            flag = 1;
            break
        end
    
        beta =  r' * ( r - r0 ) ./ rr;
        p = r + beta * p;
        iter = iter + 1; %  Count the iterations
    end
end