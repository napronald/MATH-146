function [ L2 ] = lap2d(nx,ny)
%Creates the 2-D discrete laplacian matrix (without the 1/h^2 factor)
Lx=lap1d(nx);
Ly=lap1d(ny);

Ix=speye(nx);
Iy=speye(ny);

L2=kron(Iy,Lx)+kron(Ly,Ix);

end

