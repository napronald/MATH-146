clear all; clc; close all

n=1000; % 1000x1000 matrix
[A,b] = nap.rand_sdd_matrix(n); % Randomly generated SDD matrix

x0 = zeros(length(b),1); % initialize guess 
tol= 1e-14; % tolerance
max_iter = 25000; % max number of iterations

[x,iter] = myjacobi( A, b, x0, tol, max_iter);
x_mat = A\b;

iter
relative_difference = norm(x-x_mat,inf) / norm(x,inf) %relative difference

function [x,iter] = myjacobi(A, b, x0, tol, max_iter)
iter = 0; % initialize iterations
x = x0; % define x for first iteration
y=zeros(size(x)); % initialize y
r = b - A*x0; % compute residual
while iter < max_iter
    y = x + 1./ diag(A) .* (b-A*y); % x^(k+1) = x^k + diag(A)^-1 (b-Ax^k)
    x = y; % set y (old x) as x (new x)
    r = b - A*y; % recompute residual each iteration
    if norm(r) < tol*norm(b) % stopping criteria
        break;
    end
    iter = iter + 1; % increment
    end
end