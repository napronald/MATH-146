clear all
clc
A = [1 -2 0 0 0 -1; 2 -2 3 0 0 0; 0 -2 3 1 0 0; 0 0 3 2 1 0; 0 0 0 -2 -4 -1; 0 0 0 0 1 2];
d = [-9;7;9;22;-34;17];

x = special_tri(A,length(A),d)

x_exact = [1;2;3;4;5;6];

rel_error_1 = norm(x-x_exact,1) / norm(x_exact,1)
rel_error_2 = norm(x-x_exact,2) / norm(x_exact,2)
rel_error_inf = norm(x-x_exact,inf) / norm(x_exact,inf)

function [x] = special_tri(A,n,d)
x = zeros(n,1);

% initialize (n-1 solution vector) d and t
d_prime = d(1:n-1);
t_sol(1) = A(1,n); % a_1
t_sol(n-1) = A(n-1,n); % c_n-1

% solve for v,t (n-1 solution vectors)
v = thomas(A(1:n-1,1:n-1),d_prime); % solution vector of reduced matrix
t = thomas(A(1:n-1,1:n-1),t_sol); % t = coefficents of x(n)

% using the equation from the last row + using substitution of 
% v(n-1) = x(n-1) + t(n-1)*x(n) -> x(n-1) = v(n-1) - t(n-1)*x(n)
x(n) = (d(n) - A(n,1)*v(1) - A(n,n-1)*v(n-1)) / (A(n,n) - A(n,1)*t(1) - A(n,n-1)*t(n-1));

% Using x(n) to solve for x(1) to x(n-1)
x(1:n-1) = v - x(n)*t; % x(i) + t(i)*x(n) = v(i) -> x(i) = v(i) - x(n)*t(i)
end

function [x] = thomas(A,d)
n=length(d);
if length(A)~=n || isbanded(A,1,1) ~= 1
    error('Incompatible vector lengths')
end
x=zeros(n,1);

for i=2:n
    c = A(i,i-1)/A(i-1,i-1);
    A(i,i) = A(i,i) - c*A(i-1,i);
    d(i) = d(i)-c*d(i-1);
end

x(n)=d(n)/A(n,n);
for i=n-1:-1:1
    x(i)=(d(i)-A(i,i+1)*x(i+1))/A(i,i);
end
end