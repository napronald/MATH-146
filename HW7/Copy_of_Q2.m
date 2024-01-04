clear all; clc; close all;


A = [1 2 3 4; 6 8 12 2; 1 2 3 4; 2 1 1 0];
b = [30;66;30;7];

A \ b
cond(A,2)


[U,sig,V] = svd(A)

r = 3;
sig = diag(sig);
truncated_sig = sig(1:r); % get rid of last column 

format long

% A = U*Sigma*V'

% U*Sigma*V'*x = b

z = U(:,1:r)'*b; % z = U'*b

y = z(1:end) ./ truncated_sig(1:end); % y_i = Z_i/Sigma_i

x = V(:,1:r)*y % x = V*y

residual = norm(A*x - b,2) % Compute residual

