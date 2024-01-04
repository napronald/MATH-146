clear all
clc

n=100;
Q = orth(randn(n));
D = diag(abs(randn(n, 1)));
A = Q*D*Q';
b = 10*rand(n,1);

R_t = cholesky(A);
R = transpose(R_t);

[y] = forward_sub(R_t,b);
[x] = back_sub(R,y);

x_matlab = mldivide(A,b);

validate = isequal(fix(x),fix(x_matlab));

diff=x-x_matlab;
max_rel_diff = max(abs(diff(:)))/max(abs(x_matlab))

cond_num = cond(A)
rel_matlab = norm(x-x_matlab) / norm(x_matlab)

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

function [x] = forward_sub(L,B)
n = length(B);
x = zeros(length(B),1); % Initialize solution vector

x(1) = B(1)/L(1,1); % Find first entry
for i = 2:n 
    sum = 0;
    for j = 1:i-1 
        sum = sum + L(i,j) * x(j); % Compute sum
    end
    x(i) = (B(i) - sum) / L(i,i); % Compute x entries
end
end

function [x] = back_sub(U,B)
n = length(B);
x = zeros(length(B),1); % Initialize solution vector

x(n) = B(n)/U(n,n); % find last entry
for i = n-1:-1:1 
    sum = 0;
    for j = n:-1:1 
        sum = sum + U(i,j) * x(j); % Compute sum 
    end
    x(i) = (B(i) - sum) / U(i,i); % Compute x entries
end
end

