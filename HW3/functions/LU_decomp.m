function [L, U] = LU_decomp(A)
n = size(A, 1);
L = eye(n);
U = A;
for j = 1:n-1 
    for i=j+1:n
        if U(j,j) == 0
            error('error')
        end
        L(i,j) = U(i,j) / U(j,j);
        U(i,j) = 0;
        U(i,j+1:n) = U(i,j+1:n) - L(i,j)*U(j,j+1:n);
    end
end
D = diag(U);
search = find(D==0);
if ~isempty(search)
    error('Zero found in diagonal')
end
end