clear all
clc

% A = [0.00005 1 1; 2 -1 1; 1 2 4]
% B = [2;3;3]

% form A matrix, B vector
A = 50*rand(100,100);
B = 50*rand(100,1);

% cond_num = cond(A);

x = solver(A,B); % Self built function

x_exact = mldivide(A,B); % Built in function

rel_difference = abs(x-x_exact); % Rel difference

max_norm = norm(rel_difference,inf) % Compute Max Norm

validation = isequal(fix(x),fix(x_exact)) % Validate

function x = solver(A,B)
% Using self built functions from nap.m file
[L, U] = nap.LU_decomp(A); 
y = nap.forward_sub(L,B);
x = nap.back_sub(U,y);
end
