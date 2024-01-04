clear all
clc

% Initialize Matrices
n=10;
A = 10*randn(n);
B = 10*randn(n);

% A = [1,2,3;4,5,6;7,8,0]; 
% B = [1,2,3;2,3,1;3,3,1];


x = nxn_solver(A,B) % Self-built function

x_matlab = mldivide(A,B) % built-in function

cond(A)
vali = norm(abs(x - x_matlab))


function [x] = nxn_solver(A,B) % Using function from nap.m file
    n = length(A); % Get length of Matrix A
    x = zeros(n,n); % initialize x solution vector
    for i=1:n % Loop through 1 to n using our previous functions to solve x
        [P, L, U] = nap.LUP(A);    
        PB = P*B(:,i);
        y = nap.forward_sub(L,PB);
        x(:,i) = nap.back_sub(U,y);
    end
end