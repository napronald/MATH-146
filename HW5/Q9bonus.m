clear all; clc;

% Adjust the size of our matrix
set_m = 50;
set_n = 10;

% Generating matrix with 2 linearly dependent columns
A = 10*rand(set_m,set_n);
A(:,2) = 2*A(:,1); % Making column 2 linearly dependent
A(:,4) = 2*A(:,3); % Making column 4 linearly dependent

% Values of randomly generated A,b
A;
b = 10*rand(set_m,1);

% Solving for coefficents using gs + backsub
[Q,R] = fixed_gram_schmidt(A);
y = Q'*b;
coe = nap.back_sub(R,y);

% Fix A matrix + check using backslash
col = [2 4];
A(:,col) = [];
coe_matlab = A \ b;

% Comparing our solution with matlab's backslash
cond_num = cond(A)
rel_diff = norm(coe-coe_matlab,2) / norm(coe_matlab,2) 

function [Q, R] = fixed_gram_schmidt(A) 
tol=1e-12;
index = 0;
[m,n] = size(A);

Q = zeros(m,n);
R = zeros(n,n);

for j=1:n 
    x=A(:,j); % Take the jth column
    for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j); % Compute ij entries of R, where i != j 
        x=x-R(i,j)*Q(:,i); % Compute perpendicular vector
    end
    R(j,j)=norm(x); % Define diagonal values of R
    % Modification, if norm(x) is below our tolerance (column is dependent)
    % we will discard the column/row of Q,R
    if norm(x) < tol % Here we will save the indices of dependent columns
        index = index + 1;
        discard(index) = j; % stores indices into discard array
    else
        Q(:,j) = x/R(j,j); % Computes the columns of Q
    end
    
end
R(:,discard) = []; % Eliminates the row/column of R
R(discard,:) = [];
Q(:,discard) = []; % Eliminates the column of Q
end


function [x] = forward_sub(L,B)
n = length(B);
x = zeros(length(B),1); % Initialize solution vector

x(1) = B(1)/L(1,1); % Find first entry
for i = 2:n 
    sum = 0;
    for j = 1:i-1 
        sum = sum + L(i,j) * x(j); % Compute sum
    end
    x(i) = (B(i) - sum) / L(i,i); % Compute x entries
end
end