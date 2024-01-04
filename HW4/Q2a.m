clear all
clc

% A = [2 6 0 0; 1 -1 2 0; 0 1 1 2; 0 0 -1 1];
% d = [14; -3; 9; 5];

N = 1000; % Set the size of the randomly generated tridiagonal matrix

% Using code from lab sheet 4 to generate tridiagonal matrix
a = 50*randn(1,N-1)';
b = 50*randn(1,N)';
c = 50*randn(1,N-1)';
A = diag(b) +diag(a,-1) +diag(c,1);
[Ds,k]=spdiags(A);
Z=zeros(N,N); 
B=spdiags(Ds, k,Z);

A = full(B); % Get the full A matrix

% Alternative way to generate randomly tridiagonal matrix
% A = toeplitz([50*rand() 50*rand() zeros(1,N)] , [50*rand() 50*rand() zeros(1,N)]);
d = 50*randn(N,1);

[x] = thomas_2(A,d); % Using self built function

x_matlab = mldivide(A,d); % Using built in function

cond_num = cond(A)
% vali = norm(abs(x - x_matlab)/abs(x_matlab))

rel_matlab = norm(x-x_matlab) / norm(x_matlab) % Relative Error from Matlab
% max_difference = norm(x-x_mat,2) / norm(x_mat,2) % Compute max difference


% rel_error = abs(x-x_matlab) / abs(x_matlab)

diff=x-x_matlab;
max_rel_diff = max(abs(diff(:)))/max(abs(x_matlab))
% 10^-16 * cond(A) = approximately 10^-13


validate = isequal(fix(x),fix(x_matlab)) % validate 

function [x] = thomas_2(A,d) % 
n=length(d); % Checking if input matrix is valid
if length(A)~=n || isbanded(A,1,1) ~= 1
    error('Incompatible vector lengths')
end
x=zeros(n,1);

for i=2:n %
    c = A(i,i-1)/A(i-1,i-1);
    A(i,i) = A(i,i) - c*A(i-1,i);
    d(i) = d(i)-c*d(i-1);
end

x(n)=d(n)/A(n,n);
for i=n-1:-1:1
    x(i)=(d(i)-A(i,i+1)*x(i+1))/A(i,i);
end
end