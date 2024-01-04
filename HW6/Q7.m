clear all; clc; close all

f = @(x,y) -1*exp((x-0.25).^2 + (y-0.25).^2);

% parameters
tol = 1e-8;
xmin =0;
ymin = 0;
L = 1;
N = 2^2;

% Verify solution for 5b
[u,n_it] = nap.FDJacobi_2D(f, tol,N, xmin,ymin, L );
n_it;
u = reshape(u,N*N,1);


%% Forward Difference for 2D Gauss Seidel
clc 
tic
for i=3:6 %for 2^3,2^4,..2^6
N = 2^i;
[u,n_it] = nap.FDGaussSeidel_2D(f, tol,N, xmin,ymin, L);
u = reshape(u,N*N,1);
num_points(i-2) = N;
num_iter(i-2) = n_it;
end
toc
Number_Gridpoints = num_points';
Number_Iterations = num_iter';
T = table(Number_Gridpoints,Number_Iterations)

%% Forward Difference for 2D Successive Over-Relaxation
tic
for i=3:9 % for 2^3,2^4,..2^9
N = 2^i;
[u,n_it] = nap.FDSOR_2D(f, tol,N, xmin,ymin, L );
u = reshape(u,N*N,1);
num_points(i-2) = N;
num_iter(i-2) = n_it;
end
toc
Number_Gridpoints = num_points';
Number_Iterations = num_iter';
T = table(Number_Gridpoints,Number_Iterations)