clear all; clc; close all

% u = @(x,y) -2*(x-y)*exp((y-0.25)^2 - (x-0.25)^2);

u = @(x,y) -2*exp((x-0.25)^2 + (y-0.25)^2);

% Parameters
N = (2^3)-1;
a=0;
b=1;
h = 1 / (N + 1) 

% Interior Gridpoints
xg=h*(1:N);
yg=h*(1:N); 
[xg,yg]=ndgrid(xg,yg);
bmat = u(xg,yg); 

% Formatting into solution vector form
b = reshape(bmat,N*N,1)

% Formulate A matrix
% I = speye(N,N);
% E = sparse(2:N,1:N-1,1,N,N);
% D = E+E'-2*I;
% A = kron(D,I)+kron(I,D);
% A = 1/h^2*full(A); % supposed to multiply by 1/h^2?

[ L2 ] = nap.lap2d(N,N)
A = (1/h^2)*full(L2);

u = A\b
