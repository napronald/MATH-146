%% using JACOBI Method 
clc;
% INPUTS 
n = 1000  ;             % Size of matrix                

[A,b] = rand_sdd_matrix(n)  ;         % Ensure Matrix is STRICTLY DIAGONALLY DOMINANT 
x0 = zeros(length(b),1)  ;                % Initial Guess
tol = 1e-12              ;                  % Specified tolerance 
max_iter = 20000          ;                  % Maximum # of iterations allowed to perform 

tic; 
[JACOBI_x, JACOBI_iter] = nap.myjacobi(A, b, x0, tol, max_iter) ;            % MY solution using JACOBI Method 
toc; 

tic; 
MATLAB_x = A \ b  ;                                        % MATLAB's solution using MATLAB's backslash 
toc;

JACOBI_iter                   % outputs the # of ITERATIONS PERFORMED

RelDiff_Jacobi = norm((JACOBI_x - MATLAB_x), inf) / norm(MATLAB_x, inf)         % Relative Difference w.r.t MAX-NORM 


%% using GAUSS-SEIDEL Method 
% clc;
% INPUTS 
n = 1000 ;                           % Size of matrix  

[A,b] = nap.rand_sdd_matrix(n)  ;         % Ensure Matrix is STRICTLY DIAGONALLY DOMINANT 

tol = 1e-12;                                                                                                            % Specified tolerance
max_iter = 250;                             % Maximum # of iterations allowed to perform 

tic; 
[GSEIDEL_x, GSEIDEL_iter] = nap.GaussSeidel(A, b, tol, max_iter) ;               % MY solution using GAUSS-SEIDEL Method
toc;  

tic;
MATLAB_x = A \ b  ;                                                              % MATLAB's solution using MATLAB's backslash 
toc;


GSEIDEL_iter                   % outputs the # of ITERATIONS PERFORMED


RelDiff_GSeidel = norm((GSEIDEL_x - MATLAB_x), inf) / norm(MATLAB_x, inf)           % Relative Difference w.r.t MAX-NORM 




%% using LU DECOMPOSITING w/ PIVOTING 
% clc;

n = 1000 ; 

A = randi([0 10], n, n)  ;
b = randi([0 10], n, 1)  ;

tic;
PIVOTING_x = nap.GE_PartialPivoting(A, b);                                          % MY solution using LU DECOMP. w/ PIVOTING 
toc;

tic;
MATLAB_x = A \ b;                                                              % MATLAB's solution using MATLAB's backslash 
toc;

RelDiff_Pivoting = norm((PIVOTING_x - MATLAB_x), inf) / norm(MATLAB_x, inf)          % Relative Difference w.r.t MAX-NORM 
