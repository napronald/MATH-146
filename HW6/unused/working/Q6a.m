clear all; clc; close all

n=1000;
[A,b] = nap.rand_sdd_matrix(n);

% A = [4 0 0; -2 6 1; -1 1 7];
% b = [3;9;-6];

x0 = zeros(length(b),1);
tol= 1e-12;
max_iter = 1000;

[x,iter] = myjacobi( A, b, x0, tol, max_iter);

x_mat = A\b;

iter
relative_difference = norm(x-x_mat,inf) / norm(x,inf)


function [x,iter] = myjacobi(A, b, x0, tol, max_iter)
x = x0;
iter = 1;
while iter <= max_iter
    for i = 1:length(A)
        temp = 0;
        for j = 1:length(A)
            if j ~= i
                temp = temp + A(i,j) * x0(j,1);
            end
        end
        x(i,1) = (b(i,1) - temp) / A(i,i);
    end
    iter = iter + 1;
    x0 = x;
    % check this part
    R = b- A * x0 ;
    if norm(R) < tol*norm(b)
        break
    end
end
end



% % broken code from class
% function [x,iter] = myjacobi( A, b, x0, tol, max_iter)
% iter = 1;
% x = x0;
% r = b - A.*x
% while iter <= max_iter
% for i=1:length(A)
%     y(i) = x(i) + r(i)/A(i,i);
% end
% x = y;
% r = b - A.*y;
% iter = iter + 1;
% end
% end