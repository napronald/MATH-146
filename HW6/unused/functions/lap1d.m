function [ L] = lap1d( n )
%Creates the 1-D discretized laplacian operator (without the 1/h^2 factor)
e=ones(n,1);
L=spdiags([e -2*e e], [-1 0 1], n, n);      
end

