clear all
clc
close all

x = [-2 -1 0 1 2 3];
y = [9;5;3;4;8;12];
plot(x,y,'o')


num_pts = 10;
terms = 3;

x_plot = linspace(min(x),max(x),num_pts);
coefs = lsfit(x,y,terms)


p = @(x) coefs(1) + coefs(2)*x + coefs(3)*x^2;


x = linspace(min(x),max(x),num_pts);
for i=1:length(x)
    y_approx(i) = p(x(i));
end
y_approx
hold on 
plot(x,y_approx,'r.-','MarkerSize',8)
xlabel('x')
ylabel('p(x)')
legend('Actual Points','Approximation Points')

function coefs = lsfit(t,b,n)
t = t(:);
b = b(:);
m = length(t);

A = ones(m,n);

for j=1:n-1
    A(:,j+1) = A(:,j).*t;
end

% A'b=A'A
B = A'*A;
y = A'*b;
coefs = B \ y;

% R_t = nap.cholesky(A);
% R = transpose(R_t);
% 
% [y] = nap.forward_sub(R_t,b);
% [x] = nap.back_sub(R,y)
% coefs = x;
end