clear all
clc

% A = [2 -1 0; -1 2 1; 0 -1 2];

% initialize matrix A
A = [10 5 2; 5 3 2; 2 2 3];

R_t = cholesky(A);

R = transpose(R_t);

A_verify = R_t*R;

rel_error = norm(abs(A_verify-A)/abs(A))
cond_num = cond(A)

validate = isequal(fix(A),fix(A_verify))

function [R_t] = cholesky(A)
n = size(A);
for j=1:n-1
    A(j,j) = sqrt(A(j,j));
    for i = j+1:n
        A(i,j) = A(i,j) / A(j,j);
    end
    for k = j+1:n
        for i = j+1:n
            A(i,k) = A(i,k) - A(i,j)*A(k,j);
        end
    end
end
A(n,n) = sqrt(A(n,n));
R_t = tril(A);
end