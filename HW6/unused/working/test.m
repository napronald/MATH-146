clear all; clc

n=10;

[A,b] = nap.rand_sdd_matrix(n)

D = diag(diag(A)); 
L = tril(A)- D;
U = triu(A)- D;

e= max(eig(-inv(D+L)*(U)));
if abs(e) >= 1
    disp ('Since the modulus of the largest Eigen value of iterative matrix is not less than 1') 
    disp ('this process is not convergent.')
    return
end

% x0 = zeros(length(b),1); % initialize guess 
% tol= 1e-12; % tolerance
% max_iter = 20000; % max number of iterations
% 
% [x,iter] = myjacobi( A, b, x0, tol, max_iter);
% % [x,iter] = GaussSeidel(A,b,tol,max_iter);
% 
% x_mat = A\b;
% 
% iter
% relative_difference = norm(x-x_mat,inf) / norm(x,inf) %relative difference


function [x,iter] = myjacobi(A, b, x0, tol, max_iter)
iter = 0; % initialize iterations
x = x0; % define x for first iteration
y=zeros(size(x)); % initialize y
r = b - A*x0; % compute residual
while iter < max_iter
    y = x + 1 ./ diag(A) .* (b-A*y); % x^(k+1) = x^k + diag(A)^-1 (b-Ax^k)
    x = y; % set y (old x) as x (new x)
    r = b - A*y; % recompute residual each iteration
    if norm(r) < tol*norm(b) % stopping criteria
        break;
    end
    iter = iter + 1; % increment
    end
end


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


% function [A,b] = rand_sdd_matrix(n)
% b = rand(n,1); % randomly generated b
% 
% A = 10*rand(n); % initial matrix  
% new_mat = A; % copy original matrix
% 
% for i=1:n
%     new_mat(i,i) = 0; % set diagonal entries = 0
% end
% 
% for k=1:n
%     s(k) = sum(new_mat(k,:)); % take sum of each row
% end
% 
% for i=1:n
%     A(i,i) = s(i) + 10; % change diagonal values of original matrix
% end
% end


function [A,b] = rand_sdd_matrix(n)
rng('shuffle');
A=rand(n);
B=A+A';% or A*A'
C=diag(diag(B)); % extraction of the diagonal part of B
D=B-C;
m=max(sum(D)); % m is the biggest sum in columns (or in rows) of D 
E=D+(m+0.001)*eye(n);
A = E;
b = rand(n,1);
end