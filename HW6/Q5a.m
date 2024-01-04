clear all; clc; close all

f = @(x,y) -2*(x-y)*exp((y-0.25)^2 - (x-0.25)^2); % Given f
N = 2^2; % Initialize N
h = 1 / (N + 1); % Compute h

% Interior Gridpoints
xg=h*(1:N); 
yg=h*(1:N); 
[xg,yg]=ndgrid(xg,yg);
bmat = f(xg,yg); % Evaluate at f

b = reshape(bmat,N*N,1) % reformat into vector