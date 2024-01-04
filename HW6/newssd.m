
n=5
[A,b] = rand_sdd_matrix(n)

function [A,b] = rand_sdd_matrix(n)
b = rand(n,1); % randomly generated b

A = rand(n); % initial matrix  

new_mat = A;



for k=1:n
    A(k,k) = sum(A(k,:)) +1; % take sum of each row
end
end
