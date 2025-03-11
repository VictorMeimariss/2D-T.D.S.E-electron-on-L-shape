clc;
clear;

A = gallery('poisson', 5);
b = ones(25,1); % Whatever we want
x = zeros(25,1);

tol = 10^-7; % Tolerance

% Residual
r = b - A*x;
p = r;

% Norm of rbs
norm0 = norm(r);

% Max iterations
Nmax = 1000; 

% Iteration counter
k = 0;

if norm(r) > tol*norm0

    for i = 1 : Nmax 
        
        rr = r' * r; 
        r0 = r;
        Ap = A * p;
        a = rr / (p' * Ap);
    
        x = x + a * p; % x(n+1) = x(n) + a * p 
        r = r - a * Ap;
    
        if norm(r) < tol * norm0
            break
        end
    
        beta =  r' * ( r - r0 ) ./ rr;
        p = r + beta * p;
        k = k + 1; %  Count the iterations
    end
end

disp(['Iterations: ', num2str(k)]);

func_x = conjgrad(A,b,tol);

T = table(x, func_x, 'VariableNames', {'x', 'func_x'});
disp(T);