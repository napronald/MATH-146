clear all; clc; close all

n=10;
[A,b] = nap.rand_sdd_matrix(n);

tol= 1e-15;
max_iter = 250;

[x,iter] = GaussSeidel(A,b,tol,max_iter);

x_mat = A\b;

iter
relative_difference = norm(x-x_mat,inf) / norm(x_mat,inf)

function [x,iter] = GaussSeidel(A,b,tol,max_iter)
% A = D + L + U
D = diag(diag(A)); 
L = tril(A)- D;
U = triu(A)- D;

x0 = zeros(length(A),1); % initial guess for x
iter = 0; % initialize iteration
x = x0; % define x

while iter < max_iter
    [vec1] = nap.forward_sub(D+L,U*x); % -inv(D+L)*U*x
    [vec2] = nap.forward_sub(D+L,b); % inv(D+L)*b
    y = -vec1 + vec2; % -inv(D+L)*U*x + inv(D+L)*b
    r = b - A*y;
    if norm(r) < tol*norm(b)
        break;
    end
    x = y;    
    iter = iter + 1; % increment
end
end

% iterative methods
% they are faster if you dont need many digits of accuracy
% 

% pros
%

% cons
% not as accurate 

% direct methods gives you more digitis of accuracy but may take longer
% depending on matrix also uses alot of moemory