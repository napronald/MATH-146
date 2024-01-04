clear all, close all, clc

n = 3  ;                   % Size of matrix/vector 

A = rand(n, n)    
b = rand(n , 1)  ; 

tol = 1e-12  ;             % Tolerance of choice 

x0 = zeros(n, 1)  ; 

[x, iter] = Jacobi_Iteration(A, b, x0, tol)            % MY function 


MATLAB_x = A \ b                                       % MATLAB's backslash 


%% Function to solve Ax = b using the JACOBI ITERATION METHOD
function [x, iter] = Jacobi_Iteration(A, b, x0, tol)
[m, n] = size(A)   ;    
 
for M = 1:m     
x = x0  ;                      % Initialize solution vector X
D = diag(diag(A))                       % D = Diagonal Matrix w/ the diagonal extracted from Matrix A 
        
iter = 0  ;                        % Initialize the interation counter  ------> Iteration counter:  # of iterations performed 

Error = inf  ;                     % Initialize the error to be very large 

while  Error > tol                 % STOPPING CRITERIA: use RELATIVE TOLERANCE -----> stop when || x_k+1 - x_k || < tol * || x_k ||  
    dx = D \ (b - A*x) ;
    x = x + dx         ;           % Updates the Kth iteration of x =====> Outputs the APPROX. SOLUTION VECTOR x

    Error = max(abs(dx ./ x))  ;      % Computes the RELATIVE ERROR  
  
end
  iter = iter + 1     ;            % Updates the iteration counter =====> Outputs the ITERATION COUNTER 

end


function [A,b] = rand_sdd_matrix(n)
            rng('shuffle');
            A=rand(n);
            B=A+A';% or A*A'
            C=diag(diag(B)); % extraction of the diagonal part of B
            D=B-C;
            m=max(sum(D)); % m is the biggest sum in columns (or in rows) of D 
            E=D+(m+0.001)*eye(n);
            A = E;
            b = rand(n,1);
          end
end