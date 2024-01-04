clear all; clc; 

u = @(x,y) -2*(x-y)*exp((y-0.25)^2 - (x-0.25)^2);

N = 3;
h = 1 / (N + 1);

xg=h*(1:N); 
yg=h*(1:N); 

[xg,yg]=ndgrid(xg,yg);


[ L2 ] = lap2d(N,N)

A = full(L2)

