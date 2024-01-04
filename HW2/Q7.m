clear all
clc

% randomizing matrix
A = 20*rand(1000,1000);
% A = 20*rand(2000,2000);

tic
[L, U] = LU_decompostion(A); %Self built function
toc

A_1 = L * U; %checking if L*U=A

validation = isequal(fix(A),fix(A_1)) %validate

error = norm(A_1-A,inf) %compute max error

function [L, U] = LU_decompostion(A)
n = size(A, 1); % initialize + rename variables
L = eye(n);
U = A;
for j = 1:n-1 % inputting L,U components into L,U matrices
    for i=j+1:n
        if U(j,j) == 0
            error('error')
        end
        L(i,j) = U(i,j) / U(j,j);
        U(i,j) = 0;
        U(i,j+1:n) = U(i,j+1:n) - L(i,j)*U(j,j+1:n);
    end
end
end

