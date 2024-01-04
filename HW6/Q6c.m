clear all; clc; close all

n = 1000  ;             % Size of matrix                

[A, b] = nap.rand_sdd_matrix(n) ;         % Ensure Matrix is STRICTLY DIAGONALLY DOMINANT 
MATLAB_x = A \ b  ;                       % MATLAB's solution using MATLAB's backslash 

%% using JACOBI Method 

x0 = zeros(length(b),1)         ;                % Initial Guess
tol = 1e-14                      ;                % Specified tolerance 
max_iter = 25000                 ;                % Maximum # of iterations allowed to perform 

tic; 
[JACOBI_x, JACOBI_iter] = nap.myjacobi(A, b, x0, tol, max_iter) ;            % MY solution using JACOBI Method 
toc; 

JACOBI_iter                   % outputs the # of ITERATIONS PERFORMED


RelDiff_Jacobi = norm((JACOBI_x - MATLAB_x), inf) / norm(MATLAB_x, inf)         % Relative Difference w.r.t MAX-NORM 


%% using GAUSS-SEIDEL Method 

tol = 1e-14   ;                           % Specified tolerance 
max_iter = 25000   ;                      % Maximum # of iterations allowed to perform 
   
tic; 
[GSEIDEL_x, GSEIDEL_iter] = nap.GaussSeidel(A, b, tol, max_iter) ;            % MY solution using GAUSS-SEIDEL Method
toc;  

GSEIDEL_iter                   % outputs the # of ITERATIONS PERFORMED

RelDiff_GSeidel = norm((GSEIDEL_x - MATLAB_x), inf) / norm(MATLAB_x, inf)           % Relative Difference w.r.t MAX-NORM 

%% using LU DECOMPOSITING w/ PIVOTING 

tic;
PIVOTING_x = nap.GE_PartialPivoting(A, b) ;                                          % MY solution using LU DECOMP. w/ PIVOTING 
toc;

RelDiff_Pivoting = norm((PIVOTING_x - MATLAB_x), inf) / norm(MATLAB_x, inf)          % Relative Difference w.r.t MAX-NORM 
