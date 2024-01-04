clear all
clc

A = [1,2,3,1;4,5,6,2;7,8,0,4;0,1,3,1];

[P, L, U] = LU_pivot(A)

rhs = P*A;
lhs =L*U;
validation = isequal(rhs,lhs)

function [P, L, U] = LU_pivot(A)
n=length(A);
L=eye(n); 
P= eye(n);
U=A;
for j=1:n-1
    [max_num, index] = max(abs(U(j:n,j)));
    index = index+j-1; 
    if index ~= j 
        U([index, j], :) = U([j, index], :);
        P([index, j], :) = P([j, index], :);
        L([index, j], 1:j-1) = L([j, index], 1:j-1);
    end
    for i=j+1:n
        L(i,j) = U(i,j) / U(j,j);
        U(i,j) = 0;
        U(i,j+1:n) = U(i,j+1:n) - L(i,j)*U(j,j+1:n);
    end
end
end