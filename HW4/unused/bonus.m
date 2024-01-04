clear all
clc
A = [1 -2 0 0 0 -1; 2 -2 3 0 0 0; 0 -2 3 1 0 0; 0 0 3 2 1 0; 0 0 0 -2 -4 -1; 0 0 0 0 1 2];
d = [-9;7;9;22;-34;17];

x = special_tri(A,length(A),d)

function [x] = special_tri(A,n,d)
d_prime = d(1:n-1);

v(1) = A(1,n); % a_1
v(n-1) = A(n-1,n); % c_n-1

y_prime = solver(A,n-1,d_prime);
z_prime = solver(A,n-1,v);

x(n) = (d(n)-A(n,1)*y_prime(1) - A(n,n-1)*y_prime(n-1)) / (A(n,n) - A(n,1)*z_prime(1) - A(n,n-1)*z_prime(n-1));

x_prime = y_prime - x(n)*z_prime;

x(1:n-1) = x_prime;
end

function [x] = solver(A,n,d)
[L, U] = modified_lu(A);

% Ly=b
y(1) = d(1);
for i =2:n 
    y(i) = d(i) - L(i,i-1)*y(i-1);
end

% Ux=y
x(n) = y(n) / U(n,n);
for i =n-1:-1:1 
    x(i) = (y(i)-U(i,i+1)*x(i+1))/U(i,i);
end
end

function [L, U] = modified_lu(A) % slightly modified lu function from HW #2
n = size(A, 1);
L = eye(n-1); % from n to n-1
U = A(1:n-1,1:n-1); % from n,n to n-1,n-1
for j = 1:n-1 
    for i=j+1:n-1
        L(i,j) = U(i,j) / U(j,j);
        U(i,j) = 0;
        U(i,j+1:n-1) = U(i,j+1:n-1) - L(i,j)*U(j,j+1:n-1);
    end
end
end