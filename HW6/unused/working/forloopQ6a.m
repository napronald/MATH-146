clear all; clc; close all

n=1000;
[A,b] = nap.rand_sdd_matrix(n);

x0 = zeros(length(b),1);
tol= 1e-12;
max_iter = 1000;

[x,iter] = myjacobi( A, b, x0, tol, max_iter);

x_mat = A\b;

iter
relative_difference = norm(x-x_mat,inf) / norm(x,inf)

function [x,iter] = myjacobi(A, b, x0, tol, max_iter)
iter = 0;
x = x0;
y=zeros(size(x));
r = b - A*x0;
while iter < max_iter
    for i=1:length(A)
        y(i) = x(i) + (r(i) / A(i,i));
    end
    x = y;
    r = b - A*y;
    if norm(r) < tol*norm(b)
        break;
    end
    iter = iter + 1;
    end
end