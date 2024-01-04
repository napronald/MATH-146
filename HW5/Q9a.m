%% Classical Gram-Schmidt
clear all; clc;

% Formulate A,b
A = 10*randn(50,10);
b = 10*randn(50,1);
cond_num = cond(A)

[Q, R] = gram_schmidt(A); % A=QR
c = Q'*b; % c = Q'*b
[x] = nap.back_sub(R,c); % Rx=c

% Compare with matlab's backslash
x_mat = A \ b;
rel_diff = norm((x-x_mat),2) / norm(x_mat,2)


%% Modified Gram-Schmidt
clear all; clc;

% Formulate A,b
A = 10*randn(50,10);
b = 10*randn(50,1);
cond_num = cond(A)


[Q, R] = modified_gram_schmidt(A); % A=QR
c = Q'*b; % c = Q'*b
[x] = nap.back_sub(R,c); % Rx=c

% Compare with matlab's backslash
x_mat = A \ b;
rel_diff = norm((x-x_mat),2) / norm(x_mat,2)


%% functions of CGS/MGS
function [Q, R] = gram_schmidt(A)
[m,n] = size(A);
for j=1:n
    x=A(:,j); % Take the jth column
    for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j); % Compute ij entries of R, where i != j 
        x=x-R(i,j)*Q(:,i); % Compute perpendicular vector to A(:,j)
    end
    R(j,j)=norm(x); % Define diagonal values of R
    Q(:,j) = x/R(j,j); % Computes the columns of Q
end
end

function [Q, R] = modified_gram_schmidt(A)
[m,n] = size(A);
V=A;
for i=1:n
    R(i,i) = norm(V(1:m,i));  % Define diagonal values of R
    Q(1:m,i) = V(1:m,i)/R(i,i); % Computes the columns of Q
    for j=i+1:n
        R(i,j) = Q(1:m,i)' * V(1:m,j); % Compute ij entries of R, where i != j 
        V(1:m,j) = V(1:m,j) - Q(1:m,i)*R(i,j); % Compute perpendicular vector
    end
end
end