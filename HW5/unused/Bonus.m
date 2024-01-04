clear all; clc;

% singular matrix generator
% if statement pop the row/column that is dependent
n = 10;
A = rand(n-1,n);
A(end+1,:) = sum(A);
A
A = [1 2 3; 4 5 6; 7 8 9]

b = 10*randn(n,1);

[Q,R] = gram_schmidt(A);
Q
R
% c = Q'*b(1:n-1); 
% [x] = nap.back_sub(R,c)


% verify = A \ b
% 
% rel_error = norm(x-verify) / norm(verify)



function [Q, R] = gram_schmidt(A)
tol=1e-10;
[m,n] = size(A); 

for j=1:n 
    x=A(:,j); 
    for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        x=x-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(x);
    if R(j,j) < tol
        R(:,j) = [];
        R(j,:) = [];
        Q(j,:) = [];
    else
        Q(:,j) = x/R(j,j);
    end
end
% R(:,index) = [];
% R(index,:) = [];
% Q(index,:) = [];
% Q(:,index) = [];
end