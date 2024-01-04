function [P, L, U] = LUP(A)
n=length(A);
L=eye(n); 
P= eye(n);
U=A;
for j=1:n-1
    [max_num, index] = max(abs(U(j:n,j)));
    index = index+j-1;
    if index ~= j
        U([index, j], :) = U([j, index], :);
        P([index, j], :) = P([j, index], :);
        L([index, j], 1:j-1) = L([j, index], 1:j-1);
    end
    for i=j+1:n
        L(i,j) = U(i,j) / U(j,j);
        U(i,j) = 0;
        U(i,j+1:n) = U(i,j+1:n) - L(i,j)*U(j,j+1:n);
    end
end
end