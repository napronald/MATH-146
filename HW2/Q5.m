clc 
clear all

% randomizing matrix + inputing size
n = 10; 
A = 20*rand(n,n);
L = tril(A);
B = 20*rand(n,1);

[x] = forward_sub(L,B) %Self built function

v = mldivide(L,B); % Matlab Built in function

rel = abs(x-v) / abs(v);
max_norm = norm(rel,inf) % Compute max norm

validation = isequal(fix(x),fix(v)) % validate 0 = False, 1 = True

function [x] = forward_sub(L,B)
n = length(B); % Find Size 
x = zeros(length(B),1); % Initialize x

x(1) = B(1)/L(1,1); % Compute x_1
for i = 2:n % rows
    sum = 0;
    for j = 1:i-1 % columns
        sum = sum + L(i,j) * x(j); % Compute sum to solve for x
    end
    x(i) = (B(i) - sum) / L(i,i); % Assign x_i to value
end
end