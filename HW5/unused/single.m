clear all; clc;

% n =4;
% A = rand(n-1,3);
% A(end+1,:) = sum(A);
% A

% A = 4x3 to 4x2 solution 2x1

A = 10*randn(40,30);
A(:,2) = 2*A(:,1) ;
A(:,5) = 2*A(:,4);
A
b = randn(40,1)

[Q,R] = gram_schmidt(A)
size(R)

y = Q'*b;

coe = nap.back_sub(R,y)

size(A)
A(:,2) = [];
A(:,5) = [];
size(A)
%%
coe_matlab = A \ b
%%

rel_diff = norm(coe-coe_matlab,2) / norm(coe_matlab)

% AA = Q*R

% M = AA'*AA;
% y = AA'*b;

% R_t = nap.cholesky(M);
% R = transpose(R_t);

% [z] = nap.forward_sub(R_t,y)
% [coefs] = nap.back_sub(R,z)
% 
% AA(:,n) = mat
% c = Q'*b
% [x] = nap.back_sub(R,c)

% verify = A \ b
% rel_error = norm(x-verify) / norm(verify)

function [Q, R] = gram_schmidt(A)
tol=1e-12;
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);
n_discard = 0;
for j=1:n 
    x=A(:,j);
    for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        x=x-R(i,j)*Q(:,i);
    end
%     norm(x);
    R(j,j)=norm(x);
    if norm(x) < tol
        n_discard = n_discard + 1;
        discard(n_discard) = j;
%         R(:,j) = [];
%         R(j,:) = [];
%         Q(:,j) = [];
    else
        Q(:,j) = x/R(j,j);
    end
    
end
for n=1:length(discard)
    R(:,discard) = [];
    R(discard,:) = [];
    Q(:,discard) = [];
end
end