clear all
clc

% Given A matrix + B vector
A = [20:2:36; 100:2:116; -20*ones(1,9); linspace(0.001, 0.005,9);21:2:37; 23:2:39;
24:2:40 ; linspace(-0.6, 0,9); linspace(-17,20,9)]; 
B =[1380; 4980; -900; 0.165; 1425; 1515; 1560; -9; 345];

% cond_number = cond(A);

[x] = Ax_B(A,B);

x_mat = mldivide(A,B);

x_exact = [1 2 3 4 5 6 7 8 9]'; % Given exact Solution Vector

rel_exact = norm(x-x_exact,2) / norm(x_exact,2) % Relative Error from self built function
rel_matlab = norm(x_mat-x_exact,2) / norm(x_exact,2) % Relative Error from Matlab
max_difference = norm(x-x_mat,2) / norm(x_mat,2) % Compute max difference


function [x] = Ax_B(A,B)
% Using functions from nap.m file
[P, L, U] = nap.LUP(A);

B = P*B; % formulate PB as B

% Apply Ly=B --> Ux=y
y = nap.forward_sub(L,B);
x = nap.back_sub(U,y);
end