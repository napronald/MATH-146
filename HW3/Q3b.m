clear all
clc

A = [1,2,3; 4,5,6; 7,8,0];
B = [1;2;3];

[x] = Ax_B(A,B);

x_matlab = mldivide(A,B);

validation = isequal(fix(x),fix(x_matlab))

function [x] = Ax_B(A,B)
% Using functions from nap.m file
[P, L, U] = nap.LUP(A);
B = P*B; % formulate PB as B

% Apply Ly=B --> Ux=y
y = nap.forward_sub(L,B);
x = nap.back_sub(U,y);
end