clear all
clc
A = [1 -2 0 0 0 -1; 2 -2 3 0 0 0; 0 -2 3 1 0 0; 0 0 3 2 1 0; 0 0 0 -2 -4 -1; 0 0 0 0 1 2];
d = [-9 7 9 22 -34 17];

% a = [0 1 1 -1];
% b = [2 -1 1 1];
% c = [6 2 2 0];
% d = [14; -3; 9; 5];

% A = [2 6 0 1; 1 -1 2 0; 0 1 1 2; 0 0 -1 1] %just a_1
% A = [2 6 0 1; 1 -1 2 0; 0 1 1 2; 1 0 -1 1];%with a_1+c_n
% d = [14; -3; 9; 5];

% A = [2 6 0 0; 1 -1 2 0; 0 1 1 2; 0 0 -1 1];
% d = [14; -3; 9; 5];

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

x = cyclic(n,a,b,c,d)


function [x] = cyclic(n,a,b,c,d)
v = zeros(1,n-1);
x = zeros(1,n);

x_til = zeros(1,n-1);
d_til = d(1:n-1);
y_til = zeros(1,n-1);
z_til = zeros(1,n-1);

v(1) = a(1);

v(n-1) = c(n-1);

y_til = solver(n-1,a,b,c,d_til);

z_til = solver(n-1,a,b,c,v);

x(n) = (d(n)-c(n)*y_til(1) - a(n)*y_til(n-1)) / (b(n) - c(n)*z_til(1) - a(n)*z_til(n-1));

x_til = y_til - x(n)*z_til;

x(1:n-1) = x_til;
end

function [x] = solver(n,a,b,c,d)
y = zeros(1,n);
x = zeros(1,n);

[L, U] = lu_deco(n,a,b,c);

y(1) = d(1);
for i =2:n
    y(i) = d(i) - L(i,i-1)*y(i-1);
end

x(n) = y(n) / U(n,n);

for i =n-1:-1:1
    x(i) = (y(i)-U(i,i+1)*x(i+1))/U(i,i);
end
end


function [L,U] = lu_deco(n,a,b,c)
L = eye(n,n);
U = zeros(n,n);
l = zeros(1,n);
u = zeros(1,n);
u(1) = b(1);

for i = 2:n
    l(i) = a(i) / u(i-1);
    u(i) = b(i) - l(i)*c(i-1);
end

for i=1:n
    U(i,i) = u(i);
    if i < n 
        L(i+1,i) = l(i+1);
        U(i,i+1) = c(i);
    end
end
end
