clear all
clc

% Initialize Matrix A
% A = [1,2,3;4,5,6;7,8,0];
n=10;
A = 10*randn(n);

[A_inv] = inverse(A); % built in function

I = A*A_inv; % Check that A*A_inv = identity

rel_error = norm(abs(I-eye(n)/abs(eye(n))))

cond_num = cond(A)

validate = isequal(round(I),eye(n))

function [x] = inverse(A) 
n = length(A); % Get the size of n
x = zeros(n,n); %Initialize x solution matrix 
B = eye(n); % Set B=identity and solve for x
for i=1:n % Loop through 1 to n using our previous functions to solve x
    [P, L, U] = nap.LUP(A);  
    PB = P*B(:,i); 
    y = nap.forward_sub(L,PB);
    x(:,i) = nap.back_sub(U,y);
end
end