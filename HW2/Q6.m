clc 
clear all

% randomizing matrix + inputing size
A = 20*rand(10,10);
U = triu(A);
B = 20*rand(10,1);

[x] = back_sub(U,B) %Self built function

v = mldivide(U,B); % Matlab Built in function

rel = abs(x-v) / abs(v);
max_norm = norm(rel,inf) % Compute error

validation = isequal(fix(x),fix(v)) % validate 

function [x] = back_sub(U,B)
n = length(B);
x = zeros(length(B),1);

% Check if U contains 0 in diagonal
D = diag(U);
search = find(D==0);
if ~isempty(search)
    error('error')
end

x(n) = B(n)/U(n,n); % Compute x_n
for i = n-1:-1:1 % rows
    sum = 0;
    for j = n:-1:1 % columns
        sum = sum + U(i,j) * x(j); % compute sum to solve for x
    end
    x(i) = (B(i) - sum) / U(i,i); % Assign x_i to value
end
end