function [u,n_it] = FDJacobi_2D(f, tol,N, xmin,ymin, L )
%2D Jacobi method for solving 2D Poisson equation, Lapl u = f 
% with homogeneous Dirichlet boundary conditions

%Inputs: f, a function that gives right hand side
%    tolerance
%    N, number of interior gridpoints in each direction
%    xmin, ymin for domain and L, length of domain (assumed square domain)

%outputs: solution u, as NxN matrix
%         n_int, number of iterations required

u=zeros(N,N); %initialize

h=L/(N+1);
xg=h*(1:N)+xmin;
yg=h*(1:N)+ymin;
[xg,yg]=ndgrid(xg,yg);

rhs=f(xg,yg);
%make a larger matrix with 0s to deal with the boundary terms (boundary
%terms will not be used from this matrix)
rhstouse=zeros(N+2,N+2);
rhstouse(2:N+1,2:N+1)=rhs;

reldiff2norm=1; %initialize for stopping
n_it=0;

while reldiff2norm >tol
    %make a larger matrix with 0s to deal with the boundary terms:
    utouse=zeros(N+2,N+2);
    utouse(2:N+1,2:N+1)=u;

    v=zeros(size(utouse));
    for i=2:N+1
        for j=2:N+1

            v(i,j)=1/4*(utouse(i-1,j)+utouse(i,j-1)+utouse(i+1,j)+...
                utouse(i,j+1))-h^2/4*rhstouse(i,j);
        end
    end
    v=v(2:N+1,2:N+1);
    reldiff2norm=norm(u-v,2)/norm(u,2);
    n_it=n_it+1;
    u=v;
end
end