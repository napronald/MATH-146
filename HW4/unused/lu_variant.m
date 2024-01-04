clear all
clc
A = [1 -2 0 0 0 -1; 2 -2 3 0 0 0; 0 -2 3 1 0 0; 0 0 3 2 1 0; 0 0 0 -2 -4 -1; 0 0 0 0 1 2];
d = [-9;7;9;22;-34;17];

[x] = to_vector(A,d)

function [x] = to_vector(A,d)
n= length(A);
for i = 1:n
    b(i) = A(i,i);
end
a(1) = A(1,n);
for i = 2:n
    a(i) = A(i,i-1);
end
c(n) = A(n,1);
for i =1:n-1
    c(i) = A(i,i+1);
end
x = cyclic(A,n,a,b,c,d);
end

function [x] = cyclic(A,n,a,b,c,d)
v = zeros(1,n-1);
x = zeros(1,n);

d_prime = d(1:n-1);

v(1) = a(1);

v(n-1) = c(n-1);

y_til = solver(A,n-1,d_prime);

z_til = solver(A,n-1,v);

x(n) = (d(n)-c(n)*y_til(1) - a(n)*y_til(n-1)) / (b(n) - c(n)*z_til(1) - a(n)*z_til(n-1));

x_til = y_til - x(n)*z_til;

x(1:n-1) = x_til;
end

function [x] = solver(A,n,d)
y = zeros(1,n);
x = zeros(1,n);

[L, U] = modified_lu(A);

y(1) = d(1);
for i =2:n
    y(i) = d(i) - L(i,i-1)*y(i-1);
end

x(n) = y(n) / U(n,n);

for i =n-1:-1:1
    x(i) = (y(i)-U(i,i+1)*x(i+1))/U(i,i);
end
end

function [L, U] = modified_lu(A) % slightly modified function from HW #2
n = size(A, 1);
L = eye(n-1);
U = A(1:n-1,1:n-1);
for j = 1:n-1 
    for i=j+1:n-1
        L(i,j) = U(i,j) / U(j,j);
        U(i,j) = 0;
        U(i,j+1:n-1) = U(i,j+1:n-1) - L(i,j)*U(j,j+1:n-1);
    end
end
end