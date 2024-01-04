clear all
clc
close all

x1 = [-2 -1 0 1 2 3];
y1 = [9 5 3 4 8 12];

A1 = [1 -2 4; 1 -1 1; 1 0 0; 1 1 1; 1 2 4; 1 3 9];
n = length(A1);
x = zeros(n,1);
b = [9;5;3;4;8;12];

n = length(x1);
A1 = zeros(n,n);

for j=1:n
    for i = 1:n
        A1(j,i) = x1(j)^(i-1);
    end
end

A = A1.' * A1;
b = A1.' * b;

% c = GE_PartialPivoting_VX(A,b)
c = nap.Ax_B(A,b)

p =@(x) c(1) + c(2)*x + c(3)*x.^2;

r = b - A*c

magnitude = norm(r,2)

for i = 1:n
    y_approx(i) = p(x1(i));
end

figure
scatter(x1,y1);
plot(x1, y1, 'b.', 'MarkerSize', 15);
hold on
plot(x1, p(x1), 'LineWidth', 1);
grid on;
