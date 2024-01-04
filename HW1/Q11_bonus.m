clear all
clc
close all

% e^x = e^(a) * e^m
% x = 100
% a = floor(10*x) 
% b = (x-a)/10

x = [-100:1:100];
exact = exp(x);

for i=1:length(x)
    a = floor(1000*x(i));
    b = (x(i)-(a*(0.001)));
    new_ex1 = 0.001;
    sol_1(i) = myexp(new_ex1,12);
    new_sol_1(i) = sol_1(i)^a;
    new_ex2 = b;
    sol_2(i) = myexp(new_ex2,12);
    sol(i) = new_sol_1(i) * sol_2(i);
end

abs_error = abs(exact-sol);
rel_error = abs((exact-sol)./(exact));

subplot(2,1,1)
semilogy(x,rel_error,'r-',x,abs_error,'b-')
axis([-100 100 10^-27 1])
xlabel('x')
ylabel('Error')
titletxt = sprintf('Absolute and Relative Error for n=12');
title(titletxt)
legend('Relative Error','Absolute Error')


for i=1:length(x)
    a = floor(1000*x(i));
    b = (x(i)-(a*(0.001)));
    new_ex1 = 0.001;
    sol_1(i) = myexp(new_ex1,50);
    new_sol_1(i) = sol_1(i)^a;
    new_ex2 = b;
    sol_2(i) = myexp(new_ex2,50);
    sol(i) = new_sol_1(i) * sol_2(i);
end

abs_error = abs(exact-sol);
rel_error = abs((exact-sol)./(exact));

subplot(2,1,2)
semilogy(x,rel_error,'r-',x,abs_error,'b-')
axis([-100 100 10^-27 1])
xlabel('x')
ylabel('Error')
titletxt = sprintf('Absolute and Relative Error for n=50');
title(titletxt)
legend('Relative Error','Absolute Error')


function [approx,Nterms]=myexp(x,terms);
    if x > 0
        newsum = 1;
        term = 1;
        n = 0;
        while n ~= terms-1
            n = n+1;
            term = term*x/n;
            newsum = newsum + term;
        end
        Nterms = n;
        approx = newsum;
    else
        x = abs(x);
        newsum = 1;
        term = 1;
        n = 0;
        while n ~= terms-1
            n = n+1;
            term = term*x/n;
            newsum = newsum + term;
        end
        Nterms = n;
        approx = newsum;
        approx = 1/approx;
    end
end