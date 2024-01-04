clear all; clc; close all

% Given
u = @(x,y) exp(x*y)+x^2;
x = 1;
y = 2;
h = 2.^-(2:12);

% Actual/Approx Solution
actual = (x^2+y^2)*exp(x*y)+2;
[approx] = five_point_method(u,x,y,h);

% Relative Error
for i = 1:length(h)
    rel_error(i) = norm(approx(i)-actual) / norm(actual);
end

loglog(h,rel_error,'bo-')
hold on

x = h(1:end-1);
y = rel_error(2:end);

logx = log(x)
logy = log(y)

x1 = logx(end-1)
y1 = logy(end-1)
x2 = logx(end)
y2 = logy(end)

slope = (y2 - y1) / (x2 - x1)

xpts = [x1-1 x2-1]
ypts = [y1-1 y2-1]
plot(xpts,ypts,'or-')

function [approx] = five_point_method(u,x,y,h)
    for i = 1:length(h)  
        approx(i) = (u(x-h(i),y) + u(x+h(i),y) + u(x,y-h(i)) + u(x,y+h(i)) - 4*u(x,y)) / h(i)^2;
    end
end