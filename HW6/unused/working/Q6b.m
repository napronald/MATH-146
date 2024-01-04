clear all; clc; close all

n=10;
[A,b] = nap.rand_sdd_matrix(n);

tol= 1e-12;
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
iter= 1; % initialize iteration counter
x( : , 1 ) = x0; % define x

% D+L
% U
 
[inv_mat] = nap.forward_sub(D+L,U*x)
 
inv(D+L)
% inverse(D+L)

while iter < max_iter
    x( : ,iter + 1 ) = -inv(D+L)*(U)*x( : ,iter) + inv(D+L)*b; % formula for gauss Seidel
    error = x( :, iter+1) - x( :, iter); % ||x-x0||
    if norm(error) < tol % stop criteria for tolerance
        break
    end
    norm(error);
    iter = iter + 1; % increment
end
x = x(: ,iter); % most updated x vector
end

% function [x] = inverse(A,U) 
% n = length(A); % Get the size of n
% y = zeros(n,n); %Initialize x solution matrix 
% for i=1:n % Loop through 1 to n using our previous functions to solve x
%     Ux = U*y(i);
%     y = nap.forward_sub(A,Ux);
% end
% end

% function [x] = inverse(A) 
% D = diag(diag(A));
% L = tril(A)- D;
% U = triu(A)- D;
% 
% n = length(A); % Get the size of n
% x = zeros(n,n); %Initialize x solution matrix 
% B = eye(n); % Set B=identity and solve for x
% 
% for i=1:n % Loop through 1 to n using our previous functions to solve x
%     [P, L, U] = nap.LUP(A);  
%     PB = P*B(:,i); 
%     y = nap.forward_sub(L,PB);
%     x(:,i) = nap.back_sub(U,y);
% end
% end