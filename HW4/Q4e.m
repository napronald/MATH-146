clear all
clc
close all

% Actual Data Points
x = [-2 -1 0 1 2 3];
y = [9;5;3;4;8;12];
plot(x,y,'o')

% Parameters
terms = 3;
num_pts = 6;

% formulate A matrix based on # of terms
A = ones(length(x),terms);
for j=1:terms-1
    A(:,j+1) = A(:,j).*x';
end

% Apply transformations to use Cholesky +forward/back sub
M = A'*A;
b = A'*y;

R_t = nap.cholesky(M);
R = transpose(R_t);
[z] = nap.forward_sub(R_t,b);
[coefs] = nap.back_sub(R,z);

r = y - A*coefs
magnitude = norm(r,inf)

p = @(x) coefs(1) + coefs(2)*x + coefs(3)*x^2;

x = linspace(min(x),max(x),num_pts);
for i=1:length(x)
    y_approx(i) = p(x(i));
end

hold on 
plot(x,y_approx,'r.-','MarkerSize',8)
xlabel('x')
ylabel('p(x)')
legend('Actual Points','Approximation Points')
title('For Second Order Polynomial')
