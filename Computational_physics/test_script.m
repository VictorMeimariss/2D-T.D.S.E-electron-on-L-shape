clc;
clear;

% Test script

n =300;
A = gallery('poisson', n);
b = ones(n^2,1); % Whatever we want
tol = sqrt(eps);
max_restarts = 300;
m = 150; % number of vectors in V
[x, iter, flag] = conjugate_gradient(A, b, tol);
if flag == 1
    disp(['Converged on ', num2str(iter), ' iterations']);
else
    disp('Did not converge');
end