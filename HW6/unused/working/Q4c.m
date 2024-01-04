clear all; clc; close all

upp = @(x) x - sin(pi*x) / pi^2; % actual solution

a = 0;
b = 1;
N = 2.^(3:12);

for j=1:10
h = (b - a) / (N(j) + 1);

x = a:h:b;

% A = (1/h^2) * (-2*eye(N(j)) + diag(ones(N(j)-1,1),1) + diag(ones(N(j)-1,1),-1));
A = (-2*eye(N(j)) + diag(ones(N(j)-1,1),1) + diag(ones(N(j)-1,1),-1));

% Solve for B
f = (h^2)*sin(pi*x(2:end-1))';
% f = sin(pi*x(2:end-1))';

B(1) = f(1)- (a/h^2);
for i=2:length(f)-1
    B(i) = f(i);
end
B(length(f)) = f(length(f)) - (b/h^2);

u = A\B'; %can also use thomas here

%check this part
for k=1:length(u)
    approx(k+1) = u(k);
end
approx(1) = a;
approx(end+1) = b;
approx;

rel_err(j) = norm(x-approx,2) / norm(x,2);
end

subplot(2,1,1)
plot(x,upp(x),'bo-',approx,upp(approx),'rx-');
xlabel('x');
ylabel('f(x)')
legend('Real Function','Predicted function')

x = N(1:end-1);
y = rel_err(2:end);
logx = log(x);
logy = log(y);

x1 = logx(end-1)
y1 = logy(end-1)
x2 = logx(end)
y2 = logy(end)

xpts = [x1-1 x2-1]
ypts = [y1-1 y2-1]

% plot(xpts,ypts,'-o')

slope =  (logy(end) - logy(end-1)) / (logx(end) - logx(end-1))

subplot(2,1,2)
plot(logx,logy,'r-o')
hold on
plot(xpts,ypts,'b-o')
xlabel('N')
ylabel('Relative Error')
title('Relative Error = |f(x)-f(x^*)|/|f(x)|')
legend('Relative Error')
