clear all; clc;
% Formulate A,b
A = 10*randn(50,20);
b = 10*randn(50,1);

[Q,R] = householder(A);

c = Q'*b; % c = Q'*b
[x] = nap.back_sub(R,c); % Rx=c

% Using Matlab's backslash to verify
x_mat = A \ b;
cond_num = cond(A)
rel_diff = norm((x-x_mat),2) / norm(x_mat,2)

function [Q,R] = householder(A)
[m, n] = size(A);
R = A; 
Q = eye(m,n);

for k = 1:n
    x = R(k:m,k); % kth column of R
    % Compute Reflection Vector
    if x(1) == 0 
        x(1)= norm(x,2) + x(1); % Edge case 
    else 
        x(1)= sign(x(1))*norm(x,2) + x(1); % Usual case
    end
    % Normalize vector
    x = x/norm(x,2); 
    V(k:m,k)= x/norm(x,2);
    R(k:m, k:n) = R(k:m, k:n) - 2*x*(x'*R(k:m, k:n)); % Form Upper Triangle
end
% Compute Q using V
for j=1:n
    for i=n:-1:1 
        Q(i:m,j) = Q(i:m,j) - 2*V(i:m,i)*(V(i:m,i)'*Q(i:m,j)); 
    end
end
R=R(1:n,:); % Truncate zero rows from R
end