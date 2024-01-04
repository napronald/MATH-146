%% Determine INVERSE of Matrix A (A^1 * A = I) 
function [invA] = Inv_Matrix(A) 
  %%% Input : A    ---> n x n Matrix 
  %%% Output: A^-1 ---> n x n Matrix

 n = length(A) ;         % Length of matrix A 
  
 I = eye(n) ;            % Initialize IDENTITY Matrix

 
for j = 1:n     % Loop over COLUMNS 
 
        % STEP 1: Solve for L, U, P using Gaussian Elimination W/ PARTIAL PIVOTING
        [L, U] = LU_Decomp(A) ;     

        % STEP 2: Solve L[y] = I for y using FORWARD SUB
        y = forward_sub(L, I(:,j)) ; 

        % STEP 3: Solve U[invA] = y for INVERSE A using BACKWARD SUB
        invA(:,j) = backward_sub(U, y) ;    

end


%%% (:,j) = COLUMN vector of every element in the jth column of the Matrix %%% 

